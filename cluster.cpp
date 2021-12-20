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
#include "classic.hpp"
#include "lsh.hpp"
#include "hypercube.hpp"
#include "LSHfrechet.hpp"
#include "frechetLloyds.hpp"

using namespace std;

// ./clustering -i nasd_input.csv -o outputfileLSHFrechet -c cluster.conf -update Mean Frechet -assignment LSH
// ./clustering -i nasd_input.csv -o outputLloydsfrechet -c cluster.conf -update Mean Vector -assignment Classic

// ./clustering -i nasd_input.csv -o outputLSHV -c cluster.conf -update Mean Vector -assignment LSH
// ./clustering -i nasd_input.csv -o outputHyper -c cluster.conf -update Mean Vector -assignment Hypercube
// ./clustering -i nasd_input.csv -o outputLloydsV -c cluster.conf -update Mean Vector -assignment Classic

int main(int argc, char *argv[])
{
	int k = 4;
	int L = 5;
	int N = 1;
	double R = 10000;

	string inputFile = "./";
	string outputFile = "./";
	string configurationFile = "./" ;
	string update;
	string assignment;
	bool complete = false;
	bool silhouette = false;
	int d;
	// read all given arguments 
	for (int i = 1; i < argc; i++)
	{
		string s = argv[i];
		// check if after -i aren't at the end of argv
		if( (s == "-i") and (argc > i+1))
			inputFile.append(argv[i+1]);
		// check if after -q aren't at the end of argv
		else if( (s == "-c") and (argc > i+1))
			configurationFile.append(argv[i+1]);
		else if( (s == "-o") and (argc > i+1))
			outputFile.append(argv[i+1]);
		else if( (s == "-update") and (argc > i+1) and (argc > i+2))
			update = argv[i+2];
		else if( s == "-assignment" and (argc > i+1))
			assignment = argv[i+1];
		else if( s == "-complete")
			complete = true;
		else if( s == "-silhouette")
			silhouette = true;
		// clear string to be sure that is empty in the next iteration
		s.clear();
	}
	if (inputFile == "./")
	{
		cout << "Error" << endl;
		exit(0);
	}
	else if(outputFile == "./")
	{
		cout << "Give input  output File" << endl;
		string temp;
		cin >> temp;
		outputFile.append(temp);
	}
	else if(configurationFile == "./")
	{
		cout << "Error" << endl;
		exit(0);
	}
	int K;
	int number_of_vector_hash_tables = 3; // L
	int number_of_vector_hash_functions = 4; //k
	int max_number_M_hypercube = 10;	// M
	int number_of_hypercube_dimensions = 3;	// k
	int number_of_probes = 2;	// probes

	bool clusterFound = false;
	// open configuration file with ifstream
	ifstream confile(configurationFile);
	if (confile.is_open())
	{	
		// read line of file to string
		string line;
		while (getline(confile, line))
		{
			string temp;
			bool found = false;
			// lines contains word then numbers so we separate them
			for (int i = 0; i < line.length(); i++)
			{
				// save words of the line to the string temp
				if(line[i] != ' ' and found == false)
					temp.push_back(line[i]);
				// after the word we have the number
				else if(line[i] != ' ' and found == true)
				{
					// rest of line is the number/ argument
					string number = line.substr(i);
					// cout << number << endl;
					if(temp == "number_of_clusters:")
					{
						K = stoi(number);
						clusterFound = true;
					}
					// check if after -c aren't at the end of argv
					else if(temp == "number_of_vector_hash_tables:")
						number_of_vector_hash_tables = stoi(number);
					else if(temp == "number_of_vector_hash_functions")
						number_of_vector_hash_functions = stoi(number);
					else if(temp == "max_number_M_hypercube:")
						max_number_M_hypercube = stoi(number);
					else if( temp == "number_of_hypercube_dimensions:")
						number_of_hypercube_dimensions = stoi(number);
					else if( temp == "number_of_probes:")
						number_of_probes = stoi(number);
					temp.clear();
					found = false;
				}
				else	// we found the word now we have to find the number	
					found = true;
			}
		}
		confile.close();
	}
	if(clusterFound == false)
	{
		cout << "Error in configurationFile" << endl;
		exit(-3);		
	}
	int ii = 0;
	vector<dVector> dataset;
	// open inputFile with ifstream
	ifstream file(inputFile);
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
			// get dismension of vector
			if (ii == 1)
				d = g.size();
			dVector ptrdVector = dVector(id, g);
			dataset.push_back(ptrdVector);
			g.clear();
		}
		file.close();
	}

	int n = dataset.size();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_distribution(1, n);
	std::unordered_set<int> centroids;
	int l = uniform_distribution(gen);

	// insert first random centroid for calculations
	centroids.insert(l);
	for(int t = 1; t <= (K+1); t++)
	{
		double *P = new double[n - t + 1];
		int *nonCendroid = new int[n - t + 1];
		P[0] = 0;
		// find max D(i) for all non centroids
		unsigned long long maxD = -1;
		int i = 1;
		for(int j = 0; j < n; j++) 
		{
			// checking if point is centroid 
			// if not save the index and find D(i) and use it to find P(i)
			// else go to the next point
			if(centroids.find(j + 1) == centroids.end())
			{
				// j is not a centroid
				// find D(i)
				double D = INFINITY;
				unordered_set<int>::iterator it;
				for(it = centroids.begin(); it != centroids.end(); ++it)
				{
					int c = *it;
					dVector d1 = dataset.at(c);
					dVector d2 = dataset.at(i - 1);
					vector<double> v1 =  d1.getVector();
					vector<double> v2 =  d2.getVector();
					double dist = EuclideanDistance(v1, v2);
					if (dist < D) 
						D = dist;
				}
				i++;
				if (D >= maxD) {
					maxD = D;
				}
			}
		}
		// i starts from 1 for non-centriods 
		// j starts from 0 for points
		i = 1;	// itereate through start of input dataset
		for(int j = 0; j < n; j++)
		{
			// check if vector at j is a centroid
			// if not save the index
			// find D(i) and use it for P(i) 
			// else check next point
			if (centroids.find(j+1) == centroids.end()) {
				// j is not a centroid
				// Compute D(i)
				double D = INFINITY;
				unordered_set<int>::iterator it;
				for(it = centroids.begin(); it != centroids.end(); ++it)
				{
					int c = *it;
					dVector d1 = dataset.at(c);
					dVector d2 = dataset.at(i - 1);
					vector<double> v1 =  d1.getVector();
					vector<double> v2 =  d2.getVector();
					double dist = EuclideanDistance(v1, v2);
					if (dist < D) 
						D = dist;
				}
				D /= maxD;
				P[i] = P[i - 1] + D * D;
				nonCendroid[i] = j + 1;
				i++;
			}
		}
		// choose new centroid: r 
		// chosen with probability proportional to D(r)^2
		std::uniform_real_distribution<float> floatDistribution(0, P[n - t]);
		// random uniformly distribution x [0, P(n − t)]
		// r [1, 2, ..., n − t] and P(r − 1) < x <= P(r), P(0) = 0
		float x = floatDistribution(generator);
		int left = 1;
		int right = n - t;
		int r = 0;
		// get r with binary search
		while (left <= right)
		{
 			r = (left + right) / 2;
			if (P[r-1] < x and x <= P[r])
				break;
			else if (x <= P[r-1])
				right = r - 1;
			else
				left = r + 1;
		}
		// add centroid r to set
		centroids.insert(nonCendroid[r]);
		delete[] nonCendroid;
		delete[] P;
	}
	centroids.erase(centroids.find(1));
	centroids.erase(centroids.find(l));
	// erase from set 1 and l because makes no sense that 1 is always centroid 
	vector<cluster> c;
	// create clusters
	unordered_set<int>::iterator it;
	for (it = centroids.begin(); it != centroids.end(); ++it) {
		cluster ptrCluster = cluster(dataset.at(*it-1));
		c.push_back(ptrCluster);
	}
	if(update == "Vector" and assignment == "Classic")
		Lloyds(dataset, outputFile, c, complete, silhouette, d);
	else if(update == "Vector" and assignment == "LSH")
		LSHVector(dataset.size(), number_of_vector_hash_tables, number_of_vector_hash_functions, 
			dataset, outputFile, c, d, complete, silhouette);
	else if(update == "Vector" and assignment == "Hypercube")
		hypercube(max_number_M_hypercube, number_of_hypercube_dimensions,
			number_of_probes, dataset, outputFile, c, d, complete, silhouette);
	else if(update == "Frechet" and assignment == "Classic")
		frechetLloyds(dataset, outputFile, c, complete, silhouette, d);
	else if(update == "Frechet" and assignment == "LSH")
		LSHFrechet(dataset.size(), number_of_vector_hash_tables, number_of_vector_hash_functions, 
			dataset, outputFile, c, d, complete, silhouette);	
	return 0;
}

double EuclideanDistance(std::vector<double>& v1, std::vector<double>& v2)
{
	double dist = 0 ;
	for(int i = 0 ; i < v1.size(); i++)
	{
		double temp = v1.at(i) - v2.at(i);
		dist = dist + temp*temp;
	}
	if(dist != 0)
		dist = sqrt(dist);
	return dist;
}