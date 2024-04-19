#include"user_funs.h"
#include"data_storage.h"

matrix ff_solve(matrix x, matrix ud1, matrix ud2) {
	matrix y = 0;
	DataStorage ds = DataStorage::getInstance();
	matrix xt = x;
	if(!isnan(ud2(0,0)))
	 xt = ud2[0] + x * ud2[1];

	for (int i = 0; i < 9; i++) {
		for (int i2 = 0; i2 < 21; i2++) {
			matrix s = napr_upl(xt, 0.05 * i2, ds.e_dot[i], ds.t[i]);
			y = y + pow(ds.ro[i][i2] - s, 2);
		}
	}

	double c = 1e10;
	if (xt(0) < 1.)
		y = y + c * (pow(1. - xt(0), 2));
	if (xt(0) > 1000.)
		y = y + c * (pow(xt(0) - 1000., 2));
	if (xt(1) < 0.)
		y = y + c * (pow(0. - xt(1), 2));
	if (xt(1) > 1.)
		y = y + c * (pow(xt(1) - 1., 2));
	if (xt(2) < 0.)
		y = y + c * (pow(0. - xt(2), 2));
	if (xt(2) > 1.)
		y = y + c * (pow(xt(2) - 1., 2));
	if (xt(3) < 1.)
		y = y + c * (pow(1. - xt(3), 2));
	if (xt(3) > 10000.)
		y = y + c * (pow(xt(3) - 10000., 2));
	if (xt(4) < 0.)
		y = y + c * (pow(0. - xt(4), 2));
	if (xt(4) > 1.)
		y = y + c * (pow(xt(4) - 1., 2));
	if (xt(5) < 1.)
		y = y + c * (pow(1. - xt(5), 2));
	if (xt(5) > 90000.)
		y = y + c * (pow(xt(5) - 90000., 2));
	if (xt(6) < 0.)
		y = y + c * (pow(0. - xt(6), 2));
	if (xt(6) > 1.)
		y = y + c * (pow(xt(6) - 1., 2));

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
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
)
{
	if (x >= t_cr)
		return a1 * e_dot - a2 * y * e_dot - a3 * pow(y, a8) * y_t_tcr;
	else
		return a1 * e_dot - a2 * y * e_dot;
}

DATA calc_val(DATA t, DATA tcr, int steps, DATA beg, DATA end, std::vector<DATA> y) {
	DATA td = t - tcr;
	if (td <= beg)return y[0];
	td -= beg;
	DATA h = (end - beg) / steps;
	int ind = int(floor(td / h));
	DATA rest = td - ind * h;
	DATA proc = rest / h;
	return y[ind] * proc + (1 - proc) * y[ind + 1];
}

void euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
), DATA* a) {

	DATA h = (end - beg) / steps;

	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h;
	std::cout << std::endl;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA Z = e_dot * exp(Q / (R * tK));
	DATA p_cr = -a[10] + a[11] * pow(Z, a[9]);
	std::cout << "PCR" << p_cr << std::endl;
	DATA l = a[0] / pow(Z, a[12]);
	DATA a1 = 1 / (b * l);
	DATA a2 = a[1] * pow(e_dot, -a[8]) * exp(-a[2] / (R * tK));
	DATA tau = 1e6 * u * b * b / 2;
	DATA a3 = a[3] * tau / D * exp(-a[4] / (R * tK));
	DATA t_cr = 100000;
	bool crit = true;

	for (int i = 1; i <= steps; i++) {
		ys[i] = ys[i - 1] + h * yprim(beg + (i - 1) * h, ys[i - 1], calc_val(beg + (i - 1) * h, t_cr, steps, beg, end, ys),
			a1, a2, a3, a[7], p_cr, t_cr);
		std::cout << ys[i] << " " << std::endl;
		if (ys[i] > p_cr && crit) {
			t_cr = beg + (i - 1) * h;
			crit = false;
		}
	}


	std::cout << "\nEuler Wartosc funkcji w punkcie " << ys[ys.size() - 1];

	std::ofstream plik("data.txt");
	for (int i = 0; i < ys.size(); i++)
		plik << ys[i] << std::endl;

	plik.close();
}