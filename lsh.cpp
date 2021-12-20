#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <bits/stdc++.h>
#include <numeric>
#include <map>
#include <time.h>
#include <math.h>

#include "lsh.hpp"
#include "cluster.hpp"
#include "classic.hpp"

using namespace std;

void LSHVector(int vSize, int L, int k, std::vector<dVector>& dataset, 
	std::string& outputFile, std::vector<cluster> clust ,int d, bool complete, bool silhouette)
{
	// create random natural number w
	// if it is between 300 and 400 we have better results
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(300, 400);
	int w = dis(gen);

	int i;
	int TableSize = vSize/4;
	vector<dVector> HashTable[L][TableSize];
	std::vector<double> V[L][k];
	double t[L][k];

	std::random_device normal_generator;
	normal_distribution<double> normal_distribution(0.0, 1.0);

	std::random_device generator;
	uniform_real_distribution<double> uni_distribution(0.0, w-1);

	// create vector with real coordinates distributed
	// according to the normal distribution
	for(int i = 0; i < L; i++)
	{	
		for(int j = 0; j < k; j++)
		{
			for(int ii = 0; ii < d; ii++)
			{
				double number = normal_distribution(normal_generator);
				V[i][j].push_back(number);
			}
			t[i][j] = uni_distribution(generator);
		}
	}

	std::uniform_int_distribution<> dis2(0, INT_MAX);

	// crete random vector rk with values from 0 to int max
	std::vector<int> r;
	for(int i = 0 ; i < k ; i++)
		r.push_back(dis2(gen));
	// insert each point to hash tables
	for(int i = 0; i < vSize; i++)
	{
		vector<int> hi[L];
		for(int j = 0; j < L; j++)
		{
			std::vector<int> h;
			for( int ii = 0; ii < k; ii++)
			{
				int hashF = hp(dataset.at(i), w, d, V[j][ii], t[j][ii]);
				h.push_back(hashF);
			}
			hi[j] = h;
		}
		for(int j = 0; j < L; j++)
		{
			int bucket = gp(r, hi[j], TableSize);
			dVector temp = dataset.at(i);
			temp.index = i;
			HashTable[j][bucket].push_back(temp);
		}
	}
	ofstream output;
	output.open(outputFile.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open())
	{
		cout << "Output Error" << endl;
		exit(-1);
	}
	clock_t begin = clock();

	// calculate min dist of centroids to use as first radii
	double minDist = INFINITY;
	for(int i = 0; i < clust.size(); i++)
	{
		dVector c1 = clust.at(i).getCentre();
		vector<double> v1 = c1.getVector();
		for(int j = i + 1; j < clust.size(); j++)
		{	
			dVector c2 = clust.at(j).getCentre();
			vector<double> v2 = c2.getVector();
			double dist = EuclideanDistance(v1, v2);
			if(dist < minDist)
				minDist = dist;
		}
	}
	double radii = minDist / 2 ;
	// flag is true when most centroids don't get new points
	bool flag = false;
	while(flag == false)
	{
		int counter = 0;
		for(int i = 0 ; i < clust.size(); i++)
		{
			dVector c = clust.at(i).getCentre();
			vector<double> v1 = c.getVector();

			vector<int> hi[L];
			for(int j = 0; j < L; j++)
			{
				std::vector<int> h;
				for( int ii = 0; ii < k; ii++)
				{
					int hashF = hp(c, w, d, V[j][ii], t[j][ii]);
					h.push_back(hashF);
				}
				hi[j] = h;
			}
			vector<int> buckets;
			for(int j = 0; j < L; j++)
			{
				int bucket = gp(r, hi[j], TableSize);
				buckets.push_back(bucket);
			}
			// flag is false when we don't add new points to centroid
			bool flag2 = false;
			for(int j = 0; j < L; j++)
			{
				int bucket = buckets.at(j);
				for(vector<dVector>::iterator itr = (HashTable[j][bucket]).begin();
					itr != (HashTable[j][bucket]).end(); itr++ )
				{
					vector<double> v2 = (itr)->getVector();
					double dist = EuclideanDistance(v1, v2);
					int index = itr->index;
					if(dist < radii && dataset.at(index).isMark() == false)
					{
						dataset.at(index).mark();
						flag2 = true;
						clust.at(i).addToCluster(*itr, d);
					}
				}

			}
			// counts in how many centroids points we added
			if (flag2 == true)
				counter++; 
		}
		radii = radii*2;
		// if in half the centrois to points were added stop
		if(counter < (k/2) )
			flag = true;
	}
	for(int i = 0; i < dataset.size(); i++)
	{
		if(dataset.at(i).isMark() == true)
			continue;
		double minC = INFINITY;
		int minI;
		for(int c = 0 ; c < clust.size(); c++)
		{
			dVector c1 = clust.at(c).getCentre();
			vector<double> v1 = c1.getVector();
			vector<double> v2 = dataset.at(i).getVector();
			double dist = EuclideanDistance(v1, v2);
			if(dist < minC)
			{
				minC = dist;
				minI = c;
			}
		}
		dataset.at(i).mark();
		clust.at(minI).addToCluster(dataset.at(i), d);
	}

	clock_t end = clock();
	double tClustering = (double)(end - begin) / CLOCKS_PER_SEC;
	output<< "Algorithm: LSH Mean Vector" << endl << endl;
	for(int i = 0; i < clust.size(); i++)
	{
		int size = clust.at(i).getClusterSize();
		output << "CLUSTER-" << i+1 << " Size: " << size << "\t";
		dVector temp = clust.at(i).getCentre();
		vector<double> d= temp.getVector();
		for(int j = 0; j < d.size(); j++)
			output << d.at(j) << " ";
		output << endl << endl;
	}
	output << "CLustering Time " << tClustering << endl;
	if(complete == true)
	{
		for(int i = 0; i < clust.size(); i++)
		{
			output << "CLUSTER-" << i+1 << "\t";
			vector<dVector> d = clust.at(i).getDataset();
			for(int j = 0; j < d.size(); j++)
			{
				string v = d.at(j).getId();
				output << v << " ";
			}
			output << endl << endl;
		}
	}
	if(silhouette == true)
	{
		// calculate silhouette takes about eternity
		for(int i = 0; i < clust.size(); i++)
		{
			output << "CLUSTER-" << i+1 << "\t";
			vector<dVector> d = clust.at(i).getDataset();
			double Si = 0;
			for(int j = 0; j < d.size(); j++)
			{
				vector<double> vectI = d.at(j).getVector();	//d2
				// find average distance of the point of the cluster
				// with the other points inside cluster
				double avgA= 0;
				for(int ii = 0; ii < d.size(); ii++)
				{
					if(ii == j)
						continue;
					vector <double> d3 = d.at(ii).getVector();
					double dist = EuclideanDistance(vectI, d3);
					avgA += (dist/d.size());
				}
				// find nearest centre of other cluster
				double minC = INFINITY;
				int minI;
				for(int c = 0 ; c < clust.size(); c++)
				{
					if(c == i)
						continue;
					dVector c1 = clust.at(c).getCentre();
					vector<double> v1 = c1.getVector();
					double dist = EuclideanDistance(v1, vectI);
					if(dist < minC)
					{
						minC = dist;
						minI = c;
					}
				}
				// calculate average distance with second best centroid
				

				vector<dVector> d2 = clust.at(minI).getDataset();

				double avgB = 0;
				for(int ii = 0; ii < d2.size(); ii++)
				{

					vector <double> d3 = d2.at(ii).getVector();
					double dist = EuclideanDistance(d3, vectI);
					avgB += (dist/d2.size());
				}

				double S;
				if(avgA < avgB)
					S = (1 - avgA/avgB)/d.size();
				else if(avgA > avgB)
					S = (avgB/avgA -1)/d.size();
				else if(avgA == avgB)
					S = 0;
				Si += S;
			}
			output << "Silhoutete\t" << Si << endl << endl;
		}
	}
	return;
}

int hp(dVector& dv, int w, int d, std::vector<double>& V, double t)
{
	vector<double> p = dv.getVector();
	// calculate p * v by using inner_product
	double pv = (double)inner_product(p.begin(), p.begin() + p.size(), V.begin(), 0.0);
	double h = (pv + t) / (double(w));
	return (int)h;
}

int gp(std::vector<int>& r, std::vector<int>& g, int TableSize)
{
	long long rh = inner_product(g.begin(), g.end() + g.size(), r.begin(), 0);
	long long M = pow(2, 52) - 5;
	int temp = euclidModuloLong(rh, M);
	int bucket = euclidModulo(temp, TableSize);
	return bucket;
}

// returns positive result
int euclidModulo(int x, int y)
{
	return ( (x % y) + y ) % y;
}

int euclidModuloLong(int x, long long y)
{
	return ( (x % y) + y) % y;
}

