#include "data_storage.h"

DataStorage* DataStorage::instance = nullptr;

DataStorage::DataStorage(){

	ifstream file("data.txt");
	double f;


	ro = vector<vector<double>>(2,vector<double>(100001));
	for (int i = 0; i < 2; i++) {
		file >> f;
		t.push_back(f);
		file >> f;
		e_dot.push_back(f);

		for(int i2 = 0; i2 < 100001; i2++){
			file >> f;
			ro[i][i2] = f;
		}
	}
}

DataStorage DataStorage::getInstance()
{
	if (DataStorage::instance == nullptr)
		DataStorage::instance = new DataStorage();
	return *DataStorage::instance;
}
