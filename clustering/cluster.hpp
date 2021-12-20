#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <map>


class dVector
{
	public:
		int index;
		dVector()
		{
			// std::cout << "Defaul dVector..." << std::endl;
		}
		dVector(std::string id, std::vector<double> v)
		{
			this->id = id;
			this->v = v;
			this->clust = false;
		}
		std::string getId(){
			return this->id;
		};
		std::vector<double> getVector(){
			return this->v;
		};
		void changeVector(std::vector<double> n)
		{
			this->v.clear();
			this->v = n;
		};
		void mark()
		{
			this->clust = true;
		};
		bool isMark()
		{
			return clust;
		}
	private:
		std::string id;
		std::vector<double> v;
		bool clust;

};

class cluster
{
	public:
		cluster(dVector &c)
		{
			this->centroid = c;
			// std::vector<double> temp = c.getVector()
			// this->vec.push_back(temp);
		};
		dVector getCentre()
		{
			return this->centroid;
		};
		void assignBuck(std::vector<int> b)
		{
			this->buckets = b;
		}
		std::vector<int> getBuckets()
		{
			return buckets;
		}
		void addToCluster(dVector &dV, int d)
		{
			this->dataset.push_back(dV);
			std::vector<double> temp;
			int n = dataset.size();
			std::vector<double> sum(d, 0.0);
			// std::vector<double> sum;
			for(int i = 0; i < dataset.size(); i++)
			{
				std::vector<double> tempV = dataset.at(i).getVector();
				for(int j = 0; j < tempV.size(); j++)
				{
					sum.at(j) += tempV.at(j)/n; 
				}
				
			}
			dVector newvec("-1", sum);
			this->centroid = newvec;

		}

		void initCluster(dVector &d)
		{
			this->dataset.push_back(d);
		}
		int getClusterSize()
		{
			return dataset.size();
		}
		std::vector<dVector>& getDataset()
		{
			return this->dataset;
		}

	private:
		dVector centroid;
		std::vector<dVector> dataset;
		std::vector<int> buckets;

};

double EuclideanDistance(std::vector<double>& v1, std::vector<double>& v2);