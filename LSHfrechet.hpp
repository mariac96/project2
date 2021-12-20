#pragma once
#include "cluster.hpp"


void LSHFrechet(int vSize, int L, int k, std::vector<dVector>& dataset, 
	std::string& outputFile, std::vector<cluster> clust ,int d, bool complete, bool silhouette);
int idFunction(std::vector<int> r, std::vector<int> g, int& TableSize);
int hp2(std::vector<double>& p, int w, int d, std::vector<double>& V, double t);