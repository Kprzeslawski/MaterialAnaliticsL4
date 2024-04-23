#include "data_storage.h"

DataStorage* DataStorage::instance = nullptr;

DataStorage::DataStorage(){

	ifstream file("data_out.csv");
	double f;
	save = false;

	ro = vector<vector<double>>(9,vector<double>(1001));
	string row, word;
	std::getline(file, row);
	std::stringstream lineStream(row);

	while (std::getline(lineStream, word, ';'))
		t.push_back(std::stod(word) + 273.);

	std::getline(file, row);
	std::stringstream lineStream2(row);

	while (std::getline(lineStream2, word, ';'))
		e_dot.push_back(std::stod(word));


	for(int i2 = 0; i2 < 1001; i2++){
		std::getline(file, row);
		std::stringstream lineStream3(row);
		int i = 0;
		while (std::getline(lineStream3, word, ';')) {
			ro[i][i2] = std::stod(word);
			i++;
		}
	}

	return;
}

DataStorage DataStorage::getInstance()
{
	if (DataStorage::instance == nullptr)
		DataStorage::instance = new DataStorage();
	return *DataStorage::instance;
}

DataStorage* DataStorage::getInstance2()
{
	if (DataStorage::instance == nullptr)
		DataStorage::instance = new DataStorage();
	return DataStorage::instance;
}