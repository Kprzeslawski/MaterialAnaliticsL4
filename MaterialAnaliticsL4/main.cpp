#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include "user_funs.h"
#include "data_storage.h"
#include "opt_alg.h"

int main() {
	DataStorage dt = DataStorage::getInstance();
	
	double epsilon = 1e-9;
	int Nmax = 15000;
	//for (int i = 0; i <= 10; i++) {
		matrix x0(11, 1);
		x0(0) = 1e-3 * 10 * m2d(rand_mat());
		x0(1) = 24990. * m2d(rand_mat()) + 10.;
		x0(2) = 1e+3 * (99 * m2d(rand_mat()) + 1);
		x0(3) = 3e+9 * m2d(rand_mat());
		x0(4) = 1e+3 * (490 * m2d(rand_mat()) + 10);
		x0(5) = 10. * m2d(rand_mat());
		x0(6) = m2d(0.);//10. * m2d(rand_mat());
		x0(7) = m2d(rand_mat());
		x0(8) = 1.e+13 * m2d(rand_mat());
		x0(9) = 1.e+12 * m2d(rand_mat());
		x0(10) = m2d(rand_mat());
		cout << x0(0) << " " << x0(1) << " " << x0(2) << " " << x0(3) << " " << x0(4) << " " << x0(5) << " " << x0(6) << " " 
			<< x0(7) << " " << x0(8) << " " << x0(9) << " " << x0(10)
			<< "\n";

		solution opt;

		opt = Powell(ff_solve, x0, epsilon, Nmax);
		//cout << x0(0) << " " << x0(1) << " " << x0(2) << " " << x0(3) << " " << x0(4) << " " << x0(5) << " " << x0(6) << " " << "\n";
		cout << opt.x(0) << " " << opt.x(1) << " " << opt.x(2) << " " << opt.x(3) << " " << opt.x(4) << " " << opt.x(5) << " " << opt.x(6) << " " 
			<< opt.x(7) << " " << opt.x(8) << " " << opt.x(9) << " " << opt.x(10)
			<< "\n";
		solution::clear_calls();

		matrix y = ff_solve(opt.x);
		cout << y() << "\n";
		//if (y > 0.1)return main();

		DataStorage::getInstance2()->save = true;
		y = ff_solve(opt.x);
		DataStorage::getInstance2()->save = false;

	//}
}

