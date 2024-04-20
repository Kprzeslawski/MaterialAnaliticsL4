#include"user_funs.h"


//DATA e_dot = 1;
//DATA t_maks = 1;
//DATA tK = 575. + 273.15;

//sta³e 
DATA R = 8.314462;

//stal
DATA b = 0.25e-9;
DATA D = 30;
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
		res = euler_2(p0,0.,1./ds.e_dot[i],1000,fun_ivm,xt,ds.e_dot[i],ds.t[i] + 273.);
			
		for (int i2 = 0; i2 < 1001; i2++) {
			y = y + pow((ds.ro[i][i2] - res[i2])/ds.ro[i][i2], 2) / 9000;
			//std::cout << ds.ro[i][i2] << " " << res[i2] << " " << y << std::endl;
		}
		if (ds.save && i == 7) {
			std::cout << "SAVING... ";
			std::ofstream plik("data_counted.txt");
			for (int i2 = 0; i2 < 1001; i2++) plik << res[i2] << std::endl;
			plik.close();
		}
	}
	
	double c = 1e10;
	if (xt(0) < 0.)
		y = y + c * (pow(xt(0), 2));
	if (xt(0) > 1e-2)
		y = y + c * (pow(xt(0) - 1e-2, 2));
	if (xt(1) < 10.)
		y = y + c * (pow(10. - xt(1), 2));
	if (xt(1) > 25000.)
		y = y + c * (pow(xt(1) - 25000., 2));
	if (xt(2) < 1e3)
		y = y + c * (pow(1e3 - xt(2), 2));
	if (xt(2) > 1e5)
		y = y + c * (pow(xt(2) - 1e5, 2));
	if (xt(3) < 0.)
		y = y + c * (pow(1. - xt(3), 2));
	if (xt(3) > 3e9)
		y = y + c * (pow(xt(3) - 3e9, 2));
	if (xt(4) < 1e4)
		y = y + c * (pow(1e4 - xt(4), 2));
	if (xt(4) > 5e5)
		y = y + c * (pow(xt(4) - 5e5, 2));
	if (xt(5) < 0.)
		y = y + c * (pow(0. - xt(5), 2));
	if (xt(5) > 10.)
		y = y + c * (pow(xt(5) - 10., 2));
	if (xt(6) < 0.)
		y = y + c * (pow(0. - xt(6), 2));
	if (xt(6) > .5)
		y = y + c * (pow(xt(6) - .5, 2));
	if (xt(7) < 0.)
		y = y + c * (pow(0. - xt(7), 2));
	if (xt(7) > 1.)
		y = y + c * (pow(xt(7) - 1., 2));
	if (xt(8) < 0.)
		y = y + c * (pow(0. - xt(8), 2));
	if (xt(8) > 1e13)
		y = y + c * (pow(xt(8) - 1e13, 2));
	if (xt(9) < 0.)
		y = y + c * (pow(0. - xt(9), 2));
	if (xt(9) > 1e12)
		y = y + c * (pow(xt(9) - 1e12, 2));
	if (xt(10) < 0.)
		y = y + c * (pow(0. - xt(10), 2));
	if (xt(10) > 1.)
		y = y + c * (pow(xt(10) - 1., 2));
	if (isnan(y()))y() = 1e200;
	//std::cout << y() << std::endl;
	return y;
}

matrix napr_upl(matrix x, double e, double e_dot, double t) {
	matrix y = matrix(1,1);
	double w = exp(-x(6) * e);
	double ex_1 = exp(x(3)/(8.314 * (t + 273)));
	double ex_2 = exp(x(5)/(8.314 * (t + 273)));
	
	y(0) = (w * x(0) * pow(e, x(1)) * ex_1 + (1 - w) * (x(4)) * ex_2) * pow(e_dot, x(2));
	return y;
}

DATA fun_ivm(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr, DATA e_dot
)
{
	if (x >= t_cr) {
		//std::cout <<"DATA: " << a3 << " " << a8 << " " << y << " " << y_t_tcr << " ";
		return a1 * e_dot - a2 * y * e_dot - a3 * pow(y, a8) * y_t_tcr;
	}
	else
		return a1 * e_dot - a2 * y * e_dot;
}

DATA calc_val(DATA t, DATA tcr, int steps, DATA beg, DATA end, std::vector<DATA> y) {
	DATA td = t - tcr;
	if (td <= beg)return 0.;
	td -= beg;
	DATA h = (end - beg) / steps;
	int ind = int(floor(td / h));
	DATA rest = td - ind * h;
	DATA proc = rest / h;
	return y[ind] * proc + (1 - proc) * y[ind + 1];
}

vector<DATA> euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr, DATA e_dot
), matrix a, DATA e_dot, DATA tK) {

	DATA h = (end - beg) / steps;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA Z = e_dot * exp(Q / (R * tK));
	DATA p_cr = -a(8) + a(9) * pow(Z, a(7));
	DATA l = a(0) / pow(Z, a(10));
	DATA a1 = 1. / (b * l);
	DATA a2 = a(1) * pow(e_dot, 0. /*-a(6)*/) * exp(-a(2) / (R * tK));
	DATA tau = 1e6 * u * b * b / 2;
	DATA a3 = a(3) * tau / D * exp(-a(4) / (R * tK));
	DATA t_cr = 100000;
	bool crit = true;
	//std::cout << "PARAMS: " << Z << " " << p_cr << " " << l << " " << a1 << " " << a2 << " " << tau << " " << a3 << std::endl;
	for (int i = 1; i <= steps; i++) {
		ys[i] = ys[i - 1] + h * yprim(beg + (i - 1) * h, ys[i - 1], calc_val(beg + (i - 1) * h, t_cr, steps, beg, end, ys),
			a1, a2, a3, a(5), p_cr, t_cr, e_dot);
		if (ys[i] > p_cr && crit) {
			t_cr = beg + (i - 1) * h;
			crit = false;
		}
	}

	//std::cout << "\nEuler Wartosc funkcji w punkcie " << ys[ys.size() - 1];
	
	return ys;
}