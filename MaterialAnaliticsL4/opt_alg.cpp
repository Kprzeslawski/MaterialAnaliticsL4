#include"opt_alg.h"

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		int i = 0;
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X0.y < X1.y)
		{
			d *= -1;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d;
				return p;
			}
		}
		solution X2;
		while (true)
		{
			++i;
			X2.x = x0 + pow(alpha, i) * d;
			X2.fit_fun(ff, ud1, ud2);
			if (X2.y >= X1.y || solution::f_calls > Nmax)
				break;
			X0 = X1;
			X1 = X2;
		}
		d > 0 ? p[0] = m2d(X0.x), p[1] = m2d(X2.x) : (p[0] = m2d(X2.x), p[1] = m2d(X0.x));
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double alfa = (sqrt(5) - 1) * 0.5;
		solution A, B, C, D;
		A.x = a;
		B.x = b;
		C.x = B.x - alfa * (B.x - A.x);
		C.fit_fun(ff, ud1, ud2);
		D.x = A.x + alfa * (B.x - A.x);
		D.fit_fun(ff, ud1, ud2);
		while (true)
		{
			if (C.y < D.y)
			{
				B = D;
				D = C;
				C.x = B.x - alfa * (B.x - A.x);
				C.fit_fun(ff, ud1, ud2);
			}
			else
			{
				A = C;
				C = D;
				D.x = A.x + alfa * (B.x - A.x);
				D.fit_fun(ff, ud1, ud2);
			}
			if (B.x - A.x < epsilon)
			{
				Xopt.x = (A.x + B.x) / 2;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.x = (A.x + B.x) / 2;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int n = get_len(x0);
		matrix D = ident_mat(n), A(n, 2);
		solution X, P, h;
		X.x = x0;
		double* ab;
		while (true)
		{
			P = X;
			for (int i = 0; i < n; ++i)
			{
				A.set_col(P.x, 0);
				A.set_col(D[i], 1);
				ab = expansion(ff, 0, 1, 1.5 , Nmax, ud1,A);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
				P.x = P.x + h.x * D[i];
			}
			if (norm(P.x - X.x) < epsilon)
			{
				Xopt = X ;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = X;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
			for (int i = 0; i < n - 1; ++i)
				D.set_col(D[i+1], i);

			D.set_col(P.x - X.x, n - 1);
			A.set_col(P.x, 0);
			A.set_col(D[n - 1], 1);
			ab = expansion(ff, 0, 1, 1.5, Nmax, ud1, A);
			h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
			X.x = P.x;
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}