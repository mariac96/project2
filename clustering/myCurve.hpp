#pragma once
#include "cluster.hpp"
#include "frechetLloyds.hpp"

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

class myCurve
{
    public:
        int index;
        myCurve()
        {
            // std::cout << "Defaul dVector..." << std::endl;
        }
        myCurve(std::string id)
        {
            this->id = id;
            this->clust = false;
        }
        myCurve(std::string id, std::vector<mypoint> v)
        {
            this->id = id;
            this->v = v;
            this->clust = false;
        }
        std::string getId(){
            return this->id;
        };
        std::vector<mypoint> getVector(){
            return this->v;
        };
        void setVector(std::vector<mypoint> v)
        {
            this->v = v;
        }
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
        std::vector<mypoint> v;
        bool clust;
};

class curveCluster
{
    public:
        curveCluster(myCurve &c)
        {
            this->centroid = c;
        };
        myCurve getCentre()
        {
            return this->centroid;
        };
        void addToCluster(myCurve dV, int d)
        {
            this->dataset.push_back(dV);
            // calculate mean of centroid and new curve
            std::vector<mypoint> v1 = dV.getVector();
            std::vector<mypoint> v2 = this->centroid.getVector();
            std::vector<mypoint> v;
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
            myCurve mc("-1", v);
            this->centroid = mc;
        }
        int getClusterSize()
        {
            return dataset.size();
        }
        std::vector<myCurve> getDataset()
        {
            return this->dataset;
        }
    private:
        myCurve centroid;
        std::vector<myCurve> dataset;
};