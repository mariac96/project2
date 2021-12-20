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
#include "frechet.hpp"

using namespace std;
// ./search -i nasd_input.csv -q nasd_query.csv -o outputcontinuous -L 2 -algorithm Frechet -metric continuous -delta 5 -o2
// ./search -i nasd_input.csv -q nasd_query.csv -o outputdiscrete -L 2 -algorithm Frechet -metric discrete -delta 5 -o2

// ./search -i nasd_input.csv -q nasd_query.csv -k 5 -o outputfileHypercube -algorithm Hypercube

// ./search -i nasd_input.csv -q nasd_query.csv -o outputLSH -algorithm LSH
int main(int argc, char *argv[])
{
	int k = 0;
	int L = 5;
	double delta;
	int M = 10, probes = 2;
	string inputFile="./",queryFile="./",outputFile="./",algorithm,metric;
	int d;
	// read all given arguments
	for (int i = 1; i < argc; i++)
	{
		string s = argv[i];
		// check if after -i aren't at the end of argv
		if( (s == "-i") and (argc > i+1))
			inputFile = argv[i+1];
		// check if after -q aren't at the end of argv
		else if( (s == "-q") and (argc > i+1))
			queryFile = argv[i+1];
		else if( (s == "-o") and (argc > i+1))
			outputFile = argv[i+1];
		else if( (s == "-k") and (argc > i+1))
			k = atoi(argv[i+1]);
		else if( (s == "-L") and (argc > i+1))
			L = atoi(argv[i+1]);
		else if( (s == "-M") and (argc > i+1))
			M = atoi(argv[i+1]);
		else if( (s == "-delta") and (argc > i+1))
			delta = atof(argv[i+1]);
		else if( (s == "-probes") and (argc > i+1))
			probes = atoi(argv[i+1]);
		else if( (s == "-algorithm") and (argc > i+1))
			algorithm = argv[i+1];
		else if( (s == "-metric") and (argc > i+1))
			metric = argv[i+1];
		// clear string to be sure that is empty in the next iteration
		s.clear();
	}
	if(k == 0 && algorithm == "LSH")
		k = 4;
	else if(k == 0 && algorithm == "Hypercube")
		k = 14;
	else if(k == 0 && algorithm == "Frechet")
		k = 2;
	if (inputFile == "./")
	{
		cout << "Give input File" << endl;
		string temp;
		cin >> temp;
		inputFile.append(temp);
	}
	else if(outputFile == "./")
	{
		cout << "Give input  output File" << endl;
		string temp;
		cin >> temp;
		outputFile.append(temp);
	}
	else if(queryFile == "./")
	{
		cout << "Give query File" << endl;
		string temp;
		cin >> temp;
		queryFile.append(temp);
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

	if(algorithm == "LSH")
		LSH(dataset.size(), L, k, queryFile, dataset, outputFile, d);
	else if(algorithm == "Hypercube")
		hyperCube(dataset.size(), M, k, probes, queryFile, dataset, outputFile, d);
	else if(algorithm == "Frechet" and metric == "discrete")
		discreteFrechet(dataset.size(), L, k, queryFile, dataset, outputFile, d, delta);
	else if(algorithm == "Frechet" and metric == "continuous")
		continuousFrechet(dataset.size(), L, k, queryFile, dataset, outputFile, d, delta);
	return 0;
}
