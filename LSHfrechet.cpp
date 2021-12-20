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
#include "frechetLloyds.hpp"
#include "myCurve.hpp"
#include "LSHfrechet.hpp"

using namespace std;

void LSHFrechet(int vSize, int L, int k, std::vector<dVector>& dataset, 
	std::string& outputFile, std::vector<cluster> clust ,int d, bool complete, bool silhouette)
{
	double delta = 5.0;
	int TableSize = vSize/4;
	vector<myCurve> HashTable[L][TableSize];
	// create random natural number w
	// if it is between 300 and 400 we have better results
	std::random_device rd1;
	std::mt19937 gen1(rd1());
	std::uniform_int_distribution<int> dis1(5, 50);
	int w = dis1(gen1);
	// int w = 200;
	std::vector<double> V[L][k];
	double ti[L][k];

	std::random_device normal_generator;
	normal_distribution<double> normal_distribution1(0.0, 1.0);

	std::random_device generator;
	uniform_real_distribution<double> uni_distribution1(0.0, w-1);

	// create vector with real coordinates distributed
	// according to the normal distribution
	for(int i = 0 ; i < L ; i++)
	{	
		for(int j = 0 ; j < k ; j++)
		{
			for(int ii = 0 ; ii < 2*d ; ii++)
			{
				double number = normal_distribution1(normal_generator);
				V[i][j].push_back(number);
			}
			ti[i][j] = uni_distribution1(generator);
		}
	}	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, INT_MAX);

	// crete random vector rk with values from 0 to int max
	std::vector<int> r;
	for(int i = 0 ; i < 2*d ; i++)
		r.push_back(dis(gen));
	
	std::random_device generator1;
	uniform_real_distribution<double> uni_distribution(0.0, delta);
	// double t = uni_distribution(generator);

	vector<curveCluster> cl;
	for(int c = 0 ; c < clust.size(); c++)
	{
		dVector c1 = clust.at(c).getCentre();
		string id = c1.getId();
		vector<double> vi = c1.getVector();
		vector<double>::iterator itr;
		int counter = 1;
		vector<mypoint> v;
		for(itr = vi.begin(); itr < vi.end(); itr++)
		{
			mypoint p(counter++, *itr);
			// p.printPoint();
			v.push_back(p);
		}
		// cout<<endl<<endl;
		myCurve mc(id, v);
		curveCluster ptrCluster(mc);
		cl.push_back(mc);
	}
	vector<myCurve> curveDataset;
	// construct R^2 cuve and insert them to curveDataset
	for(int i = 0; i < dataset.size(); i++)
	{
		string id = dataset.at(i).getId();
		vector<double> vi = dataset.at(i).getVector();
		vector<double>::iterator itr;
		int counter = 1;
		vector<mypoint> v;
		for(itr = vi.begin(); itr < vi.end(); itr++)
		{
			mypoint p(counter++, *itr);
			// p.printPoint();
			v.push_back(p);
		}
		// cout<<endl<<endl;
		myCurve c(id, v);
		curveDataset.push_back(c);
	}
	for(int i = 0; i < curveDataset.size(); i++)
	{
		string id = curveDataset.at(i).getId();
		vector<mypoint> vi = curveDataset.at(i).getVector();
		vector<int> hi[L];
		for(int j = 0; j < L; j++)
		{	
			vector<mypoint>::iterator itr;
			vector<mypoint> v;
			double t = uni_distribution(generator1);
			for(itr = vi.begin(); itr < vi.end(); itr++)
			{
				// snapping
				double tempX = itr->getX();
				double tempY = itr->getY();
				double x = Gdelta(tempX, t, delta);
				double y = Gdelta(tempY, t, delta);
				mypoint p(x, y);
				v.push_back(p);
			}
			vector<double> tempV;
			// remove continuous duplicates
			for(int m = 0; m < (v.size()-1); m++)
			{	
				double x1 = v.at(m).getX();
				double y1 = v.at(m).getY();
				double x2 = v.at(m+1).getX();
				double y2 = v.at(m+1).getY();
				if(compareDouble(x1, x2) == true and compareDouble(y1, y2)==true){
					v.erase(v.begin() + (m+1));
					m--;
				}
			}
			// concat to vector
			for(int m = 0; m < (v.size()); m++)
			{
				double x1 = v.at(m).getX();
				double y1 = v.at(m).getY();
				tempV.push_back(x1);
				tempV.push_back(y1);
			}
			// padding
			tempV.resize(2*d, 0);
			std::vector<int> h;
			for( int ii = 0; ii < k; ii++)
			{
				int hashF = hp2(tempV, w, d, V[j][ii], ti[j][ii]);
				h.push_back(hashF);
			}
			hi[j] = h;
			tempV.clear();
		}
		for(int j = 0; j < L; j++)
		{
			int bucket = idFunction(r, hi[j], TableSize);
			myCurve temp = curveDataset.at(i);
			temp.index = i;
			HashTable[j][bucket].push_back(temp);
		}
	}
	clock_t begin = clock();
	// calculate min dist of centroids to use as first radii
	double minDist = INFINITY;
	for(int i = 0; i < cl.size(); i++)
	{
		myCurve c1 = cl.at(i).getCentre();
		for(int j = i + 1; j < cl.size(); j++)
		{	
			myCurve c2 = cl.at(j).getCentre();
			double dist = DFD(c1, c2);
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
		for(int i = 0 ; i < cl.size(); i++)
		{
			myCurve mc = cl.at(i).getCentre();;
			vector<mypoint> vi = mc.getVector();
			vector<int> hi[L];
			for(int j = 0; j < L; j++)
			{
				vector<mypoint>::iterator itr;
				vector<mypoint> v;
				double t = uni_distribution(generator1);
				for(itr = vi.begin(); itr < vi.end(); itr++)
				{
					// snapping
					double tempX = itr->getX();
					double tempY = itr->getY();
					double x = Gdelta(tempX, t, delta);
					double y = Gdelta(tempY, t, delta);
					mypoint p(x, y);
					v.push_back(p);
				}
				vector<double> tempV;
				// remove continuous duplicates
				for(int m = 0; m < (v.size()-1); m++)
				{	
					double x1 = v.at(m).getX();
					double y1 = v.at(m).getY();
					double x2 = v.at(m+1).getX();
					double y2 = v.at(m+1).getY();
					if(compareDouble(x1, x2) == true and compareDouble(y1, y2)==true){
						v.erase(v.begin() + (m+1));
						m--;
					}
				}
				// concat to vector
				for(int m = 0; m < (v.size()); m++)
				{
					double x1 = v.at(m).getX();
					double y1 = v.at(m).getY();
					tempV.push_back(x1);
					tempV.push_back(y1);
				}
				// padding
				tempV.resize(2*d, 0);
				std::vector<int> h;
				for( int ii = 0; ii < k; ii++)
				{
					int hashF = hp2(tempV, w, d, V[j][ii], ti[j][ii]);
					h.push_back(hashF);
				}
				hi[j] = h;
				tempV.clear();
			}
			vector<int> buckets;
			for(int j = 0; j < L; j++)
			{
				int bucket = idFunction(r, hi[j], TableSize);
				buckets.push_back(bucket);
			}	
			// flag is false when we don't add new points to centroid
			bool flag2 = false;
			for(int j = 0; j < L; j++)
			{
				int bucket = buckets.at(j);
				for(vector<myCurve>::iterator itr = (HashTable[j][bucket]).begin();
				itr != (HashTable[j][bucket]).end(); itr++ )
				{
					double dist = DFD(mc, *itr);
					int index = itr->index;
					if(dist < radii && curveDataset.at(index).isMark() == false)
					{
						curveDataset.at(index).mark();
						flag2 = true;
						cl.at(i).addToCluster(*itr, d);
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
	for(int i = 0; i < curveDataset.size(); i++)
	{
		if(curveDataset.at(i).isMark() == true)
			continue;
		double minC = INFINITY;
		int minI;
		for(int c = 0 ; c < cl.size(); c++)
		{
			myCurve c1 = cl.at(c).getCentre();
			double dist = DFD(c1, curveDataset.at(i));
			if(dist < minC)
			{
				minC = dist;
				minI = c;
			}
		}
		curveDataset.at(i).mark();
		cl.at(minI).addToCluster(curveDataset.at(i), d);
	}
	ofstream output;
	output.open(outputFile.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open())
	{
		cout << "Output Error" << endl;
		exit(-1);
	}
	output<< "Algorithm: LSH Mean Frechet" << endl << endl;
	for(int i = 0; i < cl.size(); i++)
	{
		int size = cl.at(i).getClusterSize();
		output << "CLUSTER-" << i+1 << " Size: " << size << "\t";
		myCurve temp = cl.at(i).getCentre();
		std::vector<mypoint> d = temp.getVector();
		for(int j = 0; j < d.size(); j++){
			double x = d.at(j).getX();
			double y = d.at(j).getY();
			output <<"(" << x << ", " << y << ") ";
		}
		output << endl<<endl;
	}
	clock_t end = clock();
	double tClustering = (double)(end - begin) / CLOCKS_PER_SEC;
	output << "CLustering Time " << tClustering << endl;

	if(complete == true)
	{
		output << endl;
		for(int i = 0; i < cl.size(); i++)
		{
			output << "CLUSTER-" << i+1 << "\t";
			vector<myCurve> d = cl.at(i).getDataset();
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
		for(int i = 0; i < cl.size(); i++)
		{
			output << "CLUSTER-" << i+1 << "\t";
			vector<myCurve> d = cl.at(i).getDataset();
			double Si = 0;
			for(int j = 0; j < d.size(); j++)
			{
				myCurve vectI = d.at(j);
				// find average distance of the point of the cluster
				// with the other points inside cluster
				double avgA= 0;
				for(int ii = 0; ii < d.size(); ii++)
				{
					if(ii == j)
						continue;
					myCurve d3 = d.at(ii);
					double dist = DFD(vectI, d3);
					avgA += (dist/d.size());
				}
				// find nearest centre of other cluster
				double minC = INFINITY;
				int minI;
				for(int c = 0 ; c < cl.size(); c++)
				{
					if(c == i)
						continue;
					myCurve c1 = cl.at(c).getCentre();
					// vector<mypoint> v1 = c1.getVector();
					double dist = DFD(c1, vectI);
					if(dist < minC)
					{
						minC = dist;
						minI = c;
					}
				}
				// calculate average distance with second best centroid
				vector<myCurve> d2 = cl.at(minI).getDataset();
				double avgB = 0;
				for(int ii = 0; ii < d2.size(); ii++)
				{
					// myCurve d3 = d2.at(ii);
					double dist = DFD(d2.at(ii), vectI);
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
}

int hp2(std::vector<double>& p, int w, int d, std::vector<double>& V, double t)
{
	// calculate p * v by using inner_product
	double pv = (double)inner_product(p.begin(), p.begin() + p.size(), V.begin(), 0.0);
	double h = (pv + t) / (double(w));
	return (int)h;
}

int idFunction(std::vector<int> r, std::vector<int> g, int& TableSize)
{
	long long  sum = 0;
	long long int M = pow(2, 52) - 5;
	for(int i = 0; i < g.size(); i++)
	{
		long long temp1 = euclidModuloLong(r.at(i), M);
		long long int temp2 = g.at(i);
		temp2 = euclidModuloLong(temp2, M);
		sum += euclidModuloLong(temp1*temp2, M);
	} 
	sum = euclidModulo(sum, TableSize);
	
	return sum;
}