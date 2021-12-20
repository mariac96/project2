#pragma once
#include "cluster.hpp"


void Lloyds(std::vector<dVector>& dataset, std::string& outputFile,
 std::vector<cluster> clust, bool complete, bool silhouette, int dim);