#include"user_funs.h"

//sta³e 
DATA R = 8.314462;

//stal
DATA b = 0.25e-9;
DATA d = 30;
DATA u = 43500;
DATA Q = 312000;
DATA p0 = 1e4;

matrix ff_solve(matrix x, matrix ud1, matrix ud2) {
	matrix y = 0;
	DataStorage ds = DataStorage::getInstance();
	matrix xt = x;
	if(!isnan(ud2(0,0)))
	 xt = ud2[0] + x * ud2[1];
	vector<DATA> res;
	for (int i = 0; i < 9; i++) {
		//std::cout << "________________________ITER: "<< i << std::endl;
		res = CalcUsingEuler(DT{ xt(0),xt(1),xt(2),
			xt(3),xt(4),xt(5),xt(6),xt(7),xt(8),xt(9),
			xt(10)}, ds.e_dot[i], ds.t[i] + 273.);
			
		for (int i2 = 0; i2 < 1001; i2++) {
			y() = y() + pow(1. - res[i2] / ds.ro[i][i2], 2.) / 9000.;
			//std::cout << pow(1. - res[i2] / ds.ro[i][i2], 2.) / 9000. << std::endl;
			if (i2 < 1000) {
				DATA grad_res = (res[i2] - res[i2 + 1]);
				DATA grad_y = (ds.ro[i][i2] - ds.ro[i][i2 + 1]);
				y() = y() + pow(grad_res - grad_y, 2.) * 1e-28 ;
				//std::cout << pow(grad_res - grad_y, 2.) * 1e-28 << std::endl;
			}

			//std::cout << ds.ro[i][i2] << " " << res[i2] << " " << y << std::endl;
		}
		if (ds.save && i == 8) {
			std::cout << "SAVING... \n";
			std::ofstream plik("data_counted.txt");
			for (int i2 = 0; i2 < 1001; i2++) plik << res[i2] << std::endl;
			plik.close();
			std::cout << "OK\n";
		}
	}
	res.clear();
	double c = 1e10;
	if (xt(0) < 0.)
		y() = y() + c * (pow(1e3*0.05 - xt(0), 2));
	if (xt(0) > 1e-3 * 0.15)
		y() = y() + c * (pow(xt(0) - 1e-3 * 0.15, 2));
	if (xt(1) < 15000.)
		y() = y() + c * (pow(15000. - xt(1), 2));
	if (xt(1) > 22000.)
		y() = y() + c * (pow(xt(1) - 22000., 2));
	if (xt(2) < 1e3 * 50)
		y() = y() + c * (pow(1e3*50 - xt(2), 2));
	if (xt(2) > 1e3 * 100)
		y() = y() + c * (pow(xt(2) - 1e3 * 100, 2));
	if (xt(3) < 3e9 * 0.1)
		y() = y() + c * (pow(3e9*0.1 - xt(3), 2));
	if (xt(3) > 3e9 * 0.8)
		y() = y() + c * (pow(xt(3) - 3e9*0.8, 2));
	if (xt(4) < 1e3 * 50)
		y() = y() + c * (pow(1e3*50 - xt(4), 2));
	if (xt(4) > 1e3 * 150)
		y() = y() + c * (pow(xt(4) - 1e3*150, 2));
	if (xt(5) < 0.2)
		y() = y() + c * (pow(0.2 - xt(5), 2));
	if (xt(5) > 0.8)
		y() = y() + c * (pow(xt(5) - 0.8, 2));
	if (xt(6) < 0.05)
		y() = y() + c * (pow(0.05 - xt(6), 2));
	if (xt(6) > 0.25)
		y() = y() + c * (pow(xt(6) - 0.25, 2));
	if (xt(7) < 0.1)
		y() = y() + c * (pow(0.1 - xt(7), 2));
	if (xt(7) > 0.8)
		y() = y() + c * (pow(xt(7) - 0.8, 2));
	if (xt(8) < 0.)
		y() = y() + c * (pow(xt(8), 2));
	if (xt(8) > 0)
		y() = y() + c * (pow(xt(8), 2));
	if (xt(9) < 1e13 * 0.00001)
		y() = y() + c * (pow(1e13 * 0.00001 - xt(9), 2));
	if (xt(9) > 1e13 * 0.00009)
		y() = y() + c * (pow(xt(9) - 1e13*0.00009, 2));
	if (xt(10) < 0.01)
		y() = y() + c * (pow(0.01 - xt(10), 2));
	if (xt(10) > 0.08)
		y() = y() + c * (pow(xt(10) - 0.08, 2));
	if (isnan(y()))y() = 1.e100;
	//std::cout << y << endl;
	return y;
}

DT CalcUsingEuler(DT a, D e_dot, D t_kel) {
	return CalcUsingEuler(p0, 0., 1. / e_dot, 1000, a, e_dot, t_kel);
}

DT CalcUsingEuler(D beg_ro, D beg_t, D end_t, int steps, DT a, D e_dot, D t_kel) {

	D h = (end_t - beg_t) / steps;
	DT ro_dt;
	ro_dt.resize(steps + 1);
	ro_dt[0] = beg_ro;

	D Z = e_dot * exp(Q / (R * t_kel));
	D p_cr = -a[8] + a[9] * pow(Z, a[7]);
	D l = a[0] / pow(Z, a[10]);
	D a1 = 1. / (b * l);
	D a2 = a[1] * pow(e_dot, -a[6]) * exp(-a[2] / (R * t_kel));
	D tau = 1e6 * u * b * b / 2.;
	D a3 = a[3] * tau / d * exp(-a[4] / (R * t_kel));
	D t_cr_ind = 0;
	bool crit = true;

	for (int i = 1; i <= steps; i++) {

		D t = beg_t + (i - 1) * h;
		ro_dt[i] = ro_dt[i - 1];

		if (crit) {
			ro_dt[i] += h * (a1 * e_dot - a2 * ro_dt[i - 1] * e_dot);
		}
		else {
			ro_dt[i] += h * (a1 * e_dot - a2 * ro_dt[i - 1] * e_dot 
				- a3 * pow(ro_dt[i - 1], a[5]) * ro_dt[i - t_cr_ind]);
		}

		if (crit && ro_dt[i] >= p_cr) {
			t_cr_ind = i;
			crit = false;
		}
	}

	return ro_dt;
}