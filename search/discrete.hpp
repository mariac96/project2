#pragma once
#include "lsh.hpp"

void discreteFrechet(int vSize, int L, int k, std::string& queryFile, 
	std::vector<dVector>& dataset, std::string& outputFile, int d, double delta);

void continuousFrechet(int vSize, int L, int k, std::string queryFile, 
	std::vector<dVector>& dataset, std::string& outputFile, int d, double delta);

double Gdelta(double& x, double& t, double& delta);

bool compareDouble(double& a, double& b);

int hp2(std::vector<double>& p, int w, int d, std::vector<double>& V, double t);

unsigned long long longModulo(unsigned long long x, unsigned long long y);

int idFunction(std::vector<int> r, std::vector<int> g, int& TableSize);

class mypoint
{
	public:
		mypoint()
		{
			// std::cout << "Defaul dVector..." << std::endl;
		}
		mypoint(double x, double y)
		{
			this->x = x;
			this->y = y;
		}
		void printPoint()
		{
			std::cout <<"(" <<x << ", " << y<< ") ";
		}
		double getX()
		{
			return x;
		}
		double getY(){
			return y;
		}
	private:
		double x;
		double y;
};

class curve
{
	public:
		int curveID;
		curve()
		{
			// std::cout << "Defaul dVector..." << std::endl;
		}
		curve(std::string id, std::vector<mypoint> v)
		{
			this->id = id;
			this->v = v;
		}
		std::string getId(){
			return this->id;
		};
		std::vector<mypoint> getVector(){
			return this->v;
		};
	private:
		std::string id;
		std::vector<mypoint> v;

};

double DFD(curve& c1, curve& c2);

double EuclideanDFD(mypoint& p1, mypoint& p2);