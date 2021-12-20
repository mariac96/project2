#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <bits/stdc++.h>
#include <numeric>
#include <map>
#include <time.h>
#include <math.h>

#include "main.hpp"
#include "lsh.hpp"
#include "discrete.hpp"

using namespace std;



void discreteFrechet(int vSize, int L, int k, std::string& queryFile,
	std::vector<dVector>& dataset, std::string& outputFile, int d, double delta)
{
	k = 2;
	L = 2;
	int TableSize = vSize/4;
	vector<curve> HashTable[L][TableSize];

	// create random natural number w
	// if it is between 300 and 400 we have better results
	std::random_device rd1;
	std::mt19937 gen1(rd1());
	std::uniform_int_distribution<int> dis1(100, 500);
	// int w = dis1(gen1);
	int w = 50;
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

	vector<curve> curveDataset;
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
		curve c(id, v);
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
			tempV.resize(2*d, 1000);
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
			// int myID = idFunction(r, hi[j], TableSize);
			// curveDataset.at(i).curveID = myID;
			// int bucket = euclidModulo(myID, TableSize);
			HashTable[j][bucket].push_back(curveDataset.at(i));
		}
	}
	// read query file, convert data to R^2 curves
	// and them to curve dataset
	vector<curve> queries;
	int ii = 0;
	ifstream file(queryFile);
	if (file.is_open())
	{
		// read line of file to string
		string line;
		while (getline(file, line))
		{
			ii++;
			vector<mypoint> g;
			string temp;
			int j = 0;
			string id ;
			int count2 = 0;
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
						count2++;
						double x = atof(temp.c_str());
						mypoint p(count2, x);
						g.push_back(p);
						// p.printPoint();
					}
					temp.clear();
				}
			}
			//the last number
			count2++;
			double x = atof(temp.c_str());
			mypoint p(count2, x);
			// p.printPoint();
			g.push_back(p);
			curve c(id, g);
			queries.push_back(c);
			g.clear();
		}
		file.close();
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
	// cout<< "Queries size\t" << queries.size()<<endl;
	for(int i = 0; i < queries.size(); i++)
	{
		string id = queries.at(i).getId();
		vector<mypoint> vi = queries.at(i).getVector();
		clock_t begin = clock();
		vector<int> hi[L];
		// empty map container that will contain distance and vector ID
		map<double, string> euclDis;
		for(int j = 0; j < L; j++)
		{
			vector<mypoint>::iterator itr;
			vector<mypoint> v;
			double t = uni_distribution(generator);
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
			// cout<<v.size();
			vector<double> tempV;
			// remove duplicates
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
			tempV.resize(2*d, 1000);

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
			

			// int curveID = idFunction(r, hi[j], TableSize);
			// int bucket = euclidModulo(curveID, TableSize);

			for(vector<curve>::iterator itr = (HashTable[j][bucket]).begin();
				itr != (HashTable[j][bucket]).end(); itr++ )
			{
				// if(itr->curveID == curveID)
				// {
					double dist = DFD(queries.at(i), *itr);
					euclDis.insert(pair<double, string>(dist, (itr)->getId()));
				// }
			}

		}

		clock_t end = clock();
		double tLSH = (double)(end - begin) / CLOCKS_PER_SEC;
		totalLSH += tLSH;

		clock_t begin2 = clock();
		double minDist = INFINITY;
		string minId;
		// brute force calculate nearest neighbour sad!
		for(int j = 0; j < vSize; j++)
		{
			double dist = DFD(curveDataset.at(j), queries.at(i));
			if( dist < minDist)
			{
				minDist = dist;
				minId = curveDataset.at(j).getId();
			}
		}
		clock_t end2 = clock();
		double tTrue = (double)(end2 - begin2) / CLOCKS_PER_SEC;
		totalT += tTrue;
		output << "LSH Frechet Discrete" << endl;
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
		output << "distanceTrue: " << minDist << endl;
		double tempAF = it->first / minDist;
		if(tempAF > MAF)
			MAF = tempAF;
		// output << "tApproximate: " << tLSH << endl;
		// output << "tTrue: " << tTrue << endl;
		output << endl << endl;
	}
	output<< "tApproximateAverage: "<< totalLSH/queries.size() << endl;
	output<< "tTrueAverage: "<< totalT/queries.size() << endl;
	output<< "MAF: " << MAF << endl;
	return;
}

double Gdelta(double& x, double& t, double& delta)
{
	return floor(abs(x -t)/delta + 1.0/2.0) * delta + t;
}

bool compareDouble(double& a, double& b)
{
	// right way to compare floating point numbers
	if (abs(a - b) < 1e-9)
		return true;
	else
		return false;
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

unsigned long long longModulo(unsigned long long x, unsigned long long y)
{
	return ( (x % y) + y ) % y;
}


int hp2(std::vector<double>& p, int w, int d, std::vector<double>& V, double t)
{
	// calculate p * v by using inner_product
	double pv = (double)inner_product(p.begin(), p.begin() + p.size(), V.begin(), 0.0);
	double h = (pv + t) / (double(w));
	return (int)h;
}

double DFD(curve& c1, curve& c2)
{
	vector<mypoint> v1 = c1.getVector();
	vector<mypoint> v2 = c2.getVector();

	double** c= new double*[v1.size()];
	for(int i = 0; i < v1.size(); i++)
		c[i] = new double[v2.size()];

	for(int i = 0; i < v1.size(); i++)
	{
		for(int j = 0; j < v2.size(); j++)
		{
			if(i == 0 && j == 0)
			{
				c[i][j] = EuclideanDFD(v1.at(i), v2.at(j)) ;
			}
			else if(i == 0 and j > 0)
				c[i][j] = max(c[i][j-1], EuclideanDFD(v1.at(i), v2.at(j)));
			else if(i > 0 and j == 0)
				c[i][j] = max(c[i-1][j], EuclideanDFD(v1.at(i), v2.at(j)));
			else {
				c[i][j] = max( min( min(c[i-1][j], c[i-1][j-1]), c[i][j-1]),
					EuclideanDFD(v1.at(i), v2.at(j)));
			}
		}
	}
	double result = c[v1.size()-1][v2.size()-1];
	for(int i = 0; i < v1.size(); i++)
		delete [] c[i];
	delete [] c;
	return result;
}

double EuclideanDFD(mypoint& p1, mypoint& p2)
{
	double dist = 0;
	double temp1 = p1.getX() - p2.getX();
	double temp2 = p1.getY() - p2.getY();
	dist += temp1*temp1;
	dist += temp2*temp2;
	if(dist != 0)
		dist = sqrt(dist);
	return dist;
}
