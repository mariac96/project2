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
#include "cluster.hpp"


using namespace std;

void Lloyds(std::vector<dVector>& dataset, std::string& outputFile,
 std::vector<cluster> clust, bool complete, bool silhouette, int dim)
{
	clock_t begin = clock();
	// assign each point to the closest cluster
	for(int i = 0; i< dataset.size(); i++)
	{
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
		clust.at(minI).addToCluster(dataset.at(i), dim);
	}
	ofstream output;
	output.open(outputFile.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open())
	{
		cout << "Output Error" << endl;
		exit(-1);
	}
	output<< "Algorithm: Lloyds Mean Vector" << endl << endl;
	for(int i = 0; i < clust.size(); i++)
	{
		int size = clust.at(i).getClusterSize();
		output << "CLUSTER-" << i+1 << " Size: " << size << "\t";
		dVector temp = clust.at(i).getCentre();
		vector<double> d = temp.getVector();
		for(int j = 0; j < d.size(); j++)
			output << d.at(j) << " ";
		output << endl<<endl;
	}
	clock_t end = clock();
	double tClustering = (double)(end - begin) / CLOCKS_PER_SEC;
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