#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <bits/stdc++.h>
#include <numeric>
#include <map>
#include <time.h>
#include <math.h>
#include <cstdlib>

#include "main.hpp"
#include "lsh.hpp"
#include "frechet.hpp"
#include "discrete.hpp"

using namespace std;

void continuousFrechet(int vSize, int L, int k, std::string queryFile,
	std::vector<dVector>& dataset, std::string& outputFile, int d, double delta)
{
	L = 2;
	k = 2;
	vector<Curve> curveData;
	for(int i = 0; i < dataset.size(); i++)
	{
		vector<double> vi = dataset.at(i).getVector();
		vector<double>::iterator itr;
		int counter = 1;
		Points ps(1);
		for(itr = vi.begin(); itr < vi.end(); itr++)
		{
			Point p(1);
			p.set(0, *itr);
			ps.add(p);
		}
		Curve c(ps);
		curveData.push_back(c);
	}
	int TableSize = vSize/4;
	vector<dVector> HashTable[TableSize];
	// create random natural number w
	std::random_device rd1;
	std::mt19937 gen1(rd1());
	// std::uniform_int_distribution<int> dis1(100, 500);
	// int w = dis1(gen1);
	int w = 50;
	std::vector<double> V[k];
	double ti[k];
	std::random_device normal_generator;
	normal_distribution<double> normal_distribution1(0.0, 1.0);

	std::random_device generator;
	uniform_real_distribution<double> uni_distribution1(0.0, w-1);
	// create vector with real coordinates distributed
	// according to the normal distribution
	for(int j = 0 ; j < k ; j++)
	{
		for(int ii = 0 ; ii < 2*d ; ii++)
		{
			double number = normal_distribution1(normal_generator);
			V[j].push_back(number);
		}
		ti[j] = uni_distribution1(generator);
	}
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, INT_MAX);

	// crete random vector rk with values from 0 to int max
	std::vector<int> r;
	for(int i = 0 ; i < d ; i++)
		r.push_back(dis(gen));
	std::random_device generator1;
	uniform_real_distribution<double> uni_distribution(0.0, delta);

	// filtering the curve
	double e = 10.0;
	for(int i = 0; i < dataset.size(); i++)
	{
		string id = dataset.at(i).getId();
		vector<double> vi = dataset.at(i).getVector();
		int j = 0;
		while(j < (vi.size()-2))
		{
			double x1 = vi.at(j) - vi.at(j+1);
			double x2 = vi.at(j+1) - vi.at(j+2);
			if( (abs(x1) < e) and (abs(x2) < e))
			{
				vi.erase(vi.begin() + (j+1));
				j++;
			}
			else
				j++;
		}
		dataset.at(i).changeVector(vi);
	}
	double t = uni_distribution(generator1);
	for(int i = 0; i < dataset.size(); i++)
	{
		vector<double> vi = dataset.at(i).getVector();
		vector<int> hi;
			vector<double> v;
			for(int ii = 0; ii < vi.size(); ii++)
			{
				// snapping
				double temp = Gdelta(vi.at(ii), t, delta);
				v.push_back(temp);
			}
			vector<double> tempV;
			// minima maxima
			for(int ii = 1; ii < (v.size()-1); ii++)
			{
				double minima = min(v.at(ii-1), v.at(ii+1));
				double maxima = max(v.at(ii-1), v.at(ii+1));
				if(v.at(ii) >= minima and v.at(ii) <= maxima)
					tempV.push_back(v.at(ii));
			}
			// padding
			tempV.resize(d, 1000);
			std::vector<int> h;
			for( int ii = 0; ii < k; ii++)
			{
				int hashF = hp2(tempV, w, d, V[ii], ti[ii]);
				h.push_back(hashF);
			}
			hi = h;
			tempV.clear();
			int bucket = idFunction(r, hi, TableSize);
			dataset.at(i).index = i;
			HashTable[bucket].push_back(dataset.at(i));

	}
	vector<dVector> queries;
	int ii = 0;
	ifstream file(queryFile);
	if (file.is_open())
	{
		// read line of file to string
		string line;
		while (getline(file, line))
		{
			ii++;
			vector<double> g;
			string temp;
			int j = 0;
			string id ;
			// lines contains multiple numbers so we separate them
			for (int i = 0; i < line.length(); i++)
			{
				// save each digit of number to the string
				if(line[i] != '\t')
					temp.push_back(line[i]);
				else
				{	// reach end of the single number
					// so convert it from string
					// first number of line is id
					// and it is not a part of vector
					j++;
					if (j == 1)
						id = temp;
					else
					{
						double x = atof(temp.c_str());
						g.push_back(x);
					}
					temp.clear();
				}
			}
			//the last number
			double x = atof(temp.c_str());
			g.push_back(x);
			dVector ptrdVector = dVector(id, g);
			queries.push_back(ptrdVector);
			g.clear();
		}
		file.close();
	}
	vector<Curve> curvequeries;
	for(int i = 0; i < queries.size(); i++)
	{
		vector<double> vi = queries.at(i).getVector();
		vector<double>::iterator itr;
		int counter = 1;
		Points ps(1);
		for(itr = vi.begin(); itr < vi.end(); itr++)
		{
			Point p(1);
			p.set(0, *itr);
			ps.add(p);
		}
		Curve c(ps);
		curvequeries.push_back(c);
	}
	double totalT = 0.0;
	double totalLSH = 0.0;
	double MAF = 0.0;
	ofstream output;
	output.open(outputFile.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open())
	{
		cout << "Output Error" << endl;
		exit(-1);
	}
	for(int i = 0; i < queries.size(); i++)
	{
		string id = queries.at(i).getId();
		vector<double> vi = queries.at(i).getVector();
		// filtering
		int j = 0;
		while(j < (vi.size()-2))
		{
			double x1 = vi.at(j) - vi.at(j+1);
			double x2 = vi.at(j+1) - vi.at(j+2);
			if( (abs(x1) < e) and (abs(x2) < e))
			{
				vi.erase(vi.begin() + (j+1));
				j++;
			}
			else
				j++;
		}
		queries.at(i).changeVector(vi);
		// empty map container that will contain distance and vector ID
		map<double, string> euclDis;
		clock_t begin = clock();
		vector<int> hi;
		vector<double> v;
		for(int ii = 0; ii < vi.size(); ii++)
		{
			// snapping
			double temp = Gdelta(vi.at(ii), t, delta);
			v.push_back(temp);
		}
		vector<double> tempV;
		// minima maxima
		for(int ii = 1; ii < (v.size()-1); ii++)
		{
			double minima = min(v.at(ii-1), v.at(ii+1));
			double maxima = max(v.at(ii-1), v.at(ii+1));
			if(v.at(ii) >= minima and v.at(ii) <= maxima)
				tempV.push_back(v.at(ii));
		}
		// padding
		tempV.resize(d, 0);
		std::vector<int> h;
		for( int ii = 0; ii < k; ii++)
		{
			int hashF = hp2(tempV, w, d, V[ii], ti[ii]);
			h.push_back(hashF);
		}
		hi = h;
		tempV.clear();
		int bucket = idFunction(r, hi, TableSize);
		for(vector<dVector>::iterator itr = (HashTable[bucket]).begin();
			itr != (HashTable[bucket]).end(); itr++ )
		{
			int index = itr->index; 
			Frechet::Continuous::Distance d = 
			Frechet::Continuous::distance(curvequeries.at(i), curveData.at(index));
			
			double dist = d.value;
			euclDis.insert(pair<double, string>(dist, (itr)->getId()));
		}
		clock_t end = clock();
		double tLSH = (double)(end - begin) / CLOCKS_PER_SEC;
		totalLSH += tLSH;

		// clock_t begin2 = clock();
		// double minDist = INFINITY;
		// string minId;
		// // brute force calculate nearest neighbour sad!
		// for(int j = 0; j < vSize; j++)
		// {
		// 	Frechet::Continuous::Distance d = 
		// 	Frechet::Continuous::distance(curveData.at(j), curvequeries.at(i));
		// 	double dist = d.value;
		// 	if( dist < minDist)
		// 	{
		// 		minDist = dist;
		// 		minId = dataset.at(j).getId();
		// 	}
		// }
		// clock_t end2 = clock();
		// double tTrue = (double)(end2 - begin2) / CLOCKS_PER_SEC;
		// totalT += tTrue;
		
		output << "LSH Frechet Continuous" << endl;
		output << "Query: "<< id << endl;
		int counter = 0;
		int N = 1;
		map<double, string>::iterator it;
		for (it = euclDis.begin(); it != euclDis.end(); it++)
		{
			counter++;
			output << "Approximate Nearest neighbor: " << it->second << endl;
			if (counter == N)
				break;
		}
		output << "distanceApproximate: " << it->first << endl;
		
		// output << "distanceTrue: " << minDist << endl;
		// double tempAF = it->first / minDist;
		// if(tempAF > MAF)
		// 	MAF = tempAF;
		// output << "tTrue: " << tTrue << endl;
		
		output << "tApproximate: " << tLSH << endl;
		output << endl << endl;
	}
	output<< "tApproximateAverage: "<< totalLSH/queries.size() << endl;
	
	// output<< "tTrueAverage: "<< totalT/queries.size() << endl;
	// output<< "MAF: " << MAF << endl;
	
	return;
}