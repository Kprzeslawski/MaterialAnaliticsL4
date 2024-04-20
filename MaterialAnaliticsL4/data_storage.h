#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

class DataStorage {
protected:
	static DataStorage* instance;
	DataStorage();
public:
	vector<vector<double>> ro;
	vector<double> t;
	vector<double> e_dot;
	bool save;
	static DataStorage getInstance();
	static DataStorage* getInstance2();
};