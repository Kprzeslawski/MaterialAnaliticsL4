#pragma once

#include"ode_solver.h"
#include<cmath>

#define M_PI 3.14159265358979323846
#define DATA double

void euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
), DATA* a);
DATA fun_ivm(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
);

DATA e_dot = 1;
DATA t_maks = 1;
DATA tK = 575. + 273.15;

//stale 

DATA R = 8.314462;

//miedz 

DATA b = 0.25e-9;
DATA D = 30;
DATA u = 45000;
DATA Q = 238000;
DATA p0 = 1e4;

matrix ff_solve(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix napr_upl(matrix x, double e, double e_dot, double t);