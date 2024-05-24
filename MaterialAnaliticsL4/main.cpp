#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <time.h>
#include "user_funs.h"
#include "data_storage.h"
#include "opt_alg.h"

float randomD()
{
	return (D)(rand()) / (D)(RAND_MAX);
}

int main() {
	DataStorage dt = DataStorage::getInstance();
	srand(time(NULL));
	double epsilon = 1e-9;
	int Nmax = 10000;
	matrix x_opt(11, 1);
	          
	x_opt(0) = 5.97786e-05;
	x_opt(1) = 21986.2;
	x_opt(2) = 96006.7;
	x_opt(3) = 1.57093e+09;
	x_opt(4) = 142517;
	x_opt(5) = 0.799957;
	x_opt(6) = 0.249891;
	x_opt(7) = 0.799908;
	x_opt(8) = 0;
	x_opt(9) = 2.2154e+08;
	x_opt(10) = 0.08;
	DataStorage::getInstance2()->save = true;
	matrix y_opt = ff_solve(x_opt);
	DataStorage::getInstance2()->save = false;


	matrix y_tr = 1e5;
	matrix x0(11, 1);
	matrix y(1, 1);
	DT res;
	while (true) {
		while (true) {
			x0(0) = 1e-3 * ( 0.05 + rand_mat()() * 0.10);
			x0(1) = 7000 * rand_mat()() + 15000.;
			x0(2) = 1e+3 * (50 * rand_mat()() + 50);
			x0(3) = 3e+9 * (0.1 + rand_mat()() * 0.8);//4.44e6;
			x0(4) = 1e+3 * (50 * rand_mat()() + 100);//151e3;
			x0(5) = 0.2 + 0.6 * rand_mat()();//
			x0(6) = 0.05 + 0.2 * rand_mat()();//0.178546;// *rand_mat()();
			x0(7) = 0.1 + 0.8 * rand_mat()();
			x0(8) = 0. * rand_mat()();

			x0(9) = 1.e+13 * (0.00001 + 0.00008 * rand_mat()());
			x0(10) = 0.01 + 0.08 * rand_mat()();
			y = ff_solve(x0);
			//cout << y() << "\n";
			if (y < y_tr) {
				//cout << x0(0) << " " << x0(1) << " " << x0(2) << " " << x0(3) << " " << x0(4) << " " << x0(5) << " " << x0(6) << " "
				//	<< x0(7) << " " << x0(8) << " " << x0(9) << " " << x0(10)
				//	<< "\n";
				break;
			}
		}

		solution opt = Powell(ff_solve, x0, epsilon, Nmax);
		solution::clear_calls();
		y = ff_solve(opt.x);
		if (y() < y_opt()) {
			cout << "NEW MIN\n";
			cout << opt.x(0) << " " << opt.x(1) << " " << opt.x(2) << " " << opt.x(3) << " " << opt.x(4) << " " << opt.x(5) << " " << opt.x(6) << " "
				<< opt.x(7) << " " << opt.x(8) << " " << opt.x(9) << " " << opt.x(10)
				<< "\n";
			cout << y() << "\n";
			y_opt() = y();
			DataStorage::getInstance2()->save = true;
			y = ff_solve(opt.x);
			DataStorage::getInstance2()->save = false;
		}
	}
		          
		//x0(0) = 0.000178978;
		//x0(1) = 4051.7;
		//x0(2) = 63826;
		//x0(3) = 1.89127e+09;
		//x0(4) = 173086;
		//x0(5) = 0.253083;
		//x0(6) = 0.253083;
		//x0(7) = 0.995813;
		//x0(8) = 5.78715e+12;
		//x0(9) = 6.95796e+11;
		//x0(10) = 0.0969457;

		//x0(0) = 2.1e-3;
		//x0(1) = 176.0;
		//x0(2) = 19.5e3;
		//x0(3) = 4.44e+06;
		//x0(4) = 151.e3;
		//x0(5) = 1.;
		//x0(6) = 0.;
		//x0(7) = 0.262;
		//x0(8) = 0.;
		//x0(9) = 6.05e+9;
		//x0(10) = 0.167;
		//DataStorage::getInstance2()->save = true;
		//y = ff_solve(x0);
		//DataStorage::getInstance2()->save = false;
		//cout << x0(0) << " " << x0(1) << " " << x0(2) << " " << x0(3) << " " << x0(4) << " " << x0(5) << " " << x0(6) << " " 
		//	<< x0(7) << " " << x0(8) << " " << x0(9) << " " << x0(10)
		//	<< "\n";

		//solution opt;
		//
		//opt = Powell(ff_solve, x0, epsilon, Nmax);
		//cout << x0(0) << " " << x0(1) << " " << x0(2) << " " << x0(3) << " " << x0(4) << " " << x0(5) << " " << x0(6) << " " << "\n";
		//cout << opt.x(0) << " " << opt.x(1) << " " << opt.x(2) << " " << opt.x(3) << " " << opt.x(4) << " " << opt.x(5) << " " << opt.x(6) << " " 
		//	<< opt.x(7) << " " << opt.x(8) << " " << opt.x(9) << " " << opt.x(10)
		//	<< "\n";
		//solution::clear_calls();
		
		//y = ff_solve(opt.x);
		//cout << y() << "\n";
		////if (y > 0.1)return main();
		//
		//DataStorage::getInstance2()->save = true;
		//y = ff_solve(opt.x);
		//DataStorage::getInstance2()->save = false;

	//}
}

