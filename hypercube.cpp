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
#include "hypercube.hpp"


using namespace std;


void hypercube(int M, int k, int probes, std::vector<dVector>& dataset,
	string& outputFile, std::vector<cluster>& clust, int d, bool complete, bool silhouette)
{

	// create random natural number w
	// if it is between 300 and 400 we have better results
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(300, 400);
	int w = dis(gen);
	int cubeSize = pow(2, k);
	vector<dVector> Cube[cubeSize];
	std::vector<double> V[k];
	double t[k];
	std::random_device normal_generator;
	normal_distribution<double> normal_distribution(0.0, 1.0);
	std::random_device generator;
	uniform_real_distribution<double> uni_distribution(0.0, w-1);
	// create vector with real coordinates distributed
	// according to the normal distribution
	for(int i = 0 ; i < k ; i++)
	{
		for(int j = 0 ; j < d ; j++)
		{
			double number = normal_distribution(normal_generator);
			V[i].push_back(number);
		}
		t[i] = uni_distribution(generator);
	}
	
	map<int, int> myMap;
	// insert each point
	for(int i = 0; i < dataset.size(); i++)
	{
		std::vector<int> h;
		for( int j = 0; j < k; j++)
		{
			int hashF = hp(dataset.at(i), w, d, V[j], t[j]);
			h.push_back(hashF);
		}
		string p = f(h, myMap);
		int bucket = binaryToDecimal(p);
		dVector temp = dataset.at(i);
		temp.index = i;
		Cube[bucket].push_back(temp);
	}
	ofstream output;
	output.open(outputFile.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open())
	{
		cout << "Output Error" << endl;
		exit(-1);
	}
	clock_t begin = clock();
	// calculate and save buckets of each cluster
	for(int i = 0 ; i < clust.size(); i++)
	{
		dVector c = clust.at(i).getCentre();
		vector<double> v1 = c.getVector();
		std::vector<int> h;
		for( int j = 0; j < k; j++)
		{
			int hashF = hp(c, w, d, V[j], t[j]);
			h.push_back(hashF);
		}
		vector<int> buckets;
		string p = f(h, myMap);
		int bucket = binaryToDecimal(p);
		buckets.push_back(bucket);

		for(int counterProbes = 0; counterProbes < probes; counterProbes++)
		{
			// bset is initialized with bits of specified binary string
			// size is arbitary
			bitset<100> bset(p);

			// flip function flips  bits i.e.  1 <-> 0
			// starting from end to the start
			for(int ii = 0; ii < counterProbes; ii++)
				bset.flip(ii);
			string tmp1 = bset.to_string();
			// keep only the part of the string we need
			int l1 = p.length();
			int l2 = tmp1.length();
			string r = tmp1.substr(l2-l1);
			int mybucket = binaryToDecimal(r);
			buckets.push_back(mybucket);
		}

		clust.at(i).assignBuck(buckets);
	}

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

			std::vector<int> h;
			for( int j = 0; j < k; j++)
			{
				int hashF = hp(c, w, d, V[j], t[j]);
				h.push_back(hashF);
			}
			vector<int> buckets;
			string p = f(h, myMap);
			int bucket = binaryToDecimal(p);
			buckets.push_back(bucket);

			for(int counterProbes = 0; counterProbes < probes; counterProbes++)
			{
				// bset is initialized with bits of specified binary string
				// size is arbitary
				bitset<100> bset(p);

				// flip function flips  bits i.e.  1 <-> 0
				// starting from end to the start
				for(int ii = 0; ii < counterProbes; ii++)
					bset.flip(ii);
				string tmp1 = bset.to_string();
				// keep only the part of the string we need
				int l1 = p.length();
				int l2 = tmp1.length();
				string r = tmp1.substr(l2-l1);
				int mybucket = binaryToDecimal(r);
				buckets.push_back(mybucket);
			}

			// flag is false when we don't add new points to centroid
			bool flag2 = false;

			for(int j = 0; j < buckets.size(); j++)
			{
				int bucket = buckets.at(j);
				for(vector<dVector>::iterator itr = (Cube[bucket]).begin();
					itr != (Cube[bucket]).end(); itr++)
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

	for(int i = 0; i< dataset.size(); i++)
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

	output<< "Algorithm: Hypercube Mean Vector" << endl << endl;
	for(int i = 0; i < clust.size(); i++)
	{
		int size = clust.at(i).getClusterSize();
		output << "CLUSTER-" << i+1 << " Size: " << size << "\t";
		dVector temp = clust.at(i).getCentre();
		vector<double> d= temp.getVector();
		for(int j = 0; j < d.size(); j++)
			output << d.at(j) << " ";
		output << endl<<endl;
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

string f(std::vector<int> h, std::map<int, int> &m)
{
	string temp;
	for(int i = 0; i < h.size(); i++)
	{
		static std::random_device r;
		static std::seed_seq seeds{ r(), r(), r(), r(), r(), r() };
		// Construct a meresenne twister random number generator,
		// seeded with the values we just generated. This will give
		// us better quality random numbers than the random device
		static std::mt19937 engine(seeds);
		// Notice how every variable up to this point is static. This
		// is to ensure that they're only ever initialized once during
		// the execution of the program. We don't want to generate a
		// new random number generator at each function call as we'll
		// lose the random state.
		static std::uniform_int_distribution<int> dist(0, 1);
		// Roll the dice
		int coin = dist(engine);
		// if h.at(i) is not present then toss the coin
		// because same h must have same value
		if (m.count(h.at(i)) == 0)
		{
			int coin = dist(engine);
			// m.insert(pair<int, int> (h.at(i), coin));
			temp += to_string(coin);
		}
		else
		{
			coin = m.at(h.at(i));
			temp += to_string(coin);
		}
	}
	return temp;
}

// function to convert binary to decimal
int binaryToDecimal(string n)
{
	string num = n;
	int dec_value = 0;
	// Initializing base value to 1, i.e 2^0
	int base = 1;
	int len = num.length();
	for (int i = len - 1; i >= 0; i--) {
		if (num[i] == '1')
			dec_value += base;
		base = base * 2;
	}
	return dec_value;
}