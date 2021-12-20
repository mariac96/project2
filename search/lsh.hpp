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
		}
		std::string getId(){
			return this->id;
		};
		std::vector<double> getVector(){
			return this->v;
		};
		void changeVector(std::vector<double> &n)
		{
			this->v.clear();
			this->v = n;
		};
	private:
		std::string id;
		std::vector<double> v;

};