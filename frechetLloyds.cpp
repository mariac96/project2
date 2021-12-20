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
#include "frechetLloyds.hpp"

using namespace std;

void frechetLloyds(std::vector<dVector>& dataset, std::string& outputFile,
 std::vector<cluster> clust, bool complete, bool silhouette, int dim)
{
	clock_t begin = clock();
	double delta = 5.0;
	std::random_device generator1;
	uniform_real_distribution<double> uni_distribution(0.0, delta);

	// create 
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
			v.push_back(p);
		}
		myCurve c(id, v);
		curveDataset.push_back(c);
	}
	for(int i = 0; i < curveDataset.size(); i++)
	{
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
		cl.at(minI).addToCluster(curveDataset.at(i), dim);
	}
	ofstream output;
	output.open(outputFile.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open())
	{
		cout << "Output Error" << endl;
		exit(-1);
	}
	output << "Algorithm: Lloyds Mean Frechet" << endl << endl;
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
	// for(int i = 0; i < curveDataset.size(); i++)
	// {
	// 	myCurve c3("-1");
	// 	c3 = meanCurve(curveDataset.at(i), curveDataset.at(i+1), c3);

	// 	vector<mypoint>::iterator itr;
	// 	vector<mypoint> m = c3.getVector();
	// 	for(itr = m.begin(); itr < m.end(); itr++)
	// 	{
	// 		itr->printPoint();
	// 	}
	// 	break;
	// }
}

// void curveCluster::addToCluster(myCurve &dV, int d)
// {
// 	this->dataset.push_back(dV);
	

// 	// vector<myCurve> tempCurve = this->getDataset();
// 	// int n = tempCurve.size() + 1;
// 	// random_shuffle(tempCurve.begin(), tempCurve.end());
// 	// myCurve c2("-1");

// 	// myCurve *array = new myCurve [n];
// 	// array[0] = c2;
// 	// for(int i = 1; i <= tempCurve.size(); i++)
// 	// 	array[i] = tempCurve.at(i-1);



// 	// myCurve c3("-1");
// 	// myCurve c = postOrder(1, n, array, c3);
	

// 	// delete [] array;

// 	// this->centroid = c;
// }


myCurve postOrder(int index, int size, myCurve *array, myCurve& c3)
{
	// if(index > 0)
	// {		
		// is leaf
		// cout<<index<<endl;
		if((2 * index) > size)
			return array[index];
		else
		{
			myCurve leftCurve = postOrder(2*index, size,array, c3);
			myCurve rightCurve;
			if((2 * index + 1) > size)
				myCurve rightCurve = postOrder((2 * index + 1), size, array, c3);
			else 
				return leftCurve;
			return meanCurve(leftCurve, rightCurve, c3);
		}
	// }
}

myCurve meanCurve(myCurve& c1, myCurve& c2, myCurve& c3)
{
	vector<mypoint> v1 = c1.getVector();
	vector<mypoint> v2 = c2.getVector();
	vector<mypoint> v;
	for(int i = 0; i < v1.size(); i++)
	{
		double x1 = v1.at(i).getX();
		double x2 = v2.at(i).getX();
		double y1 = v1.at(i).getY();
		double y2 = v2.at(i).getY();

		double meanX = (x1 + x2)/2;
		double meanY = (y1 + y2)/2;
		mypoint p(meanX, meanY);
		v.push_back(p);
	}
	c3.setVector(v);
	return c3;
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



double DFD(myCurve& c1, myCurve& c2)
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