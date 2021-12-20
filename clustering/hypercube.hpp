#pragma once
#include "cluster.hpp"

void hypercube(int M, int k, int probes, std::vector<dVector>& dataset,
	std::string& outputFile, std::vector<cluster>& clust, int d, bool complete, bool silhouette);

std::string f(std::vector<int> h, std::map<int, int> &m);

int binaryToDecimal(std::string n);