#pragma once
#include "cluster.hpp"
#include "myCurve.hpp"


void frechetLloyds(std::vector<dVector>& dataset, std::string& outputFile,
 std::vector<cluster> clust, bool complete, bool silhouette, int dim);



double DFD(myCurve& c1, myCurve& c2);
double EuclideanDFD(mypoint& p1, mypoint& p2);

double Gdelta(double& x, double& t, double& delta);

bool compareDouble(double& a, double& b);

myCurve meanCurve(myCurve& c1, myCurve& c2, myCurve& c3);
myCurve postOrder(int index, int size, myCurve *array, myCurve& c3);