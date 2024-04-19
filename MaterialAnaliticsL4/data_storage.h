#pragma once
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class DataStorage {
protected:
	static DataStorage* instance;
	DataStorage();
public:
	vector<vector<double>> ro;
	vector<double> t;
	vector<double> e_dot;
	static DataStorage getInstance();
};