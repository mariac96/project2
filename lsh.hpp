#pragma once
#include "cluster.hpp"



int hp(dVector& dv, int w, int d, std::vector<double>& V, double t);

int euclidModuloLong(int x, long long y);

int euclidModulo(int x, int y);

int gp(std::vector<int>& r, std::vector<int>& g, int TableSize);


void LSHVector(int vSize, int L, int k, std::vector<dVector>& dataset,
	std::string& outputFile, std::vector<cluster> c,int d, 
	bool complete, bool silhouette);
