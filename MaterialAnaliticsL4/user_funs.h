#pragma once

#include"ode_solver.h"
#include"data_storage.h"
#include<cmath>

#define M_PI 3.14159265358979323846
#define DATA double
#define D double
#define DT std::vector<D>

//vector<DATA> euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr,
//	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr, DATA e_dot
//), matrix a, DATA e_dot, DATA tK);
//DATA fun_ivm(DATA x, DATA y, DATA y_t_tcr,
//	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr, DATA e_dot
//);


matrix ff_solve(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
//matrix napr_upl(matrix x, double e, double e_dot, double t);

DT CalcUsingEuler(DT a, D e_dot, D t_kel);
DT CalcUsingEuler(D beg_ro, D beg_t, D end_t, int steps, DT a, D e_dot, D t_kel);