#include"opt_alg.h"

#include <array>
#include <vector>

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

solution fib(matrix(*ff)(matrix, matrix, matrix), const double a, const double b, double epsilon, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = b - a;
        size_t n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
        std::vector<long int> F{ 1, 1 };
        F.resize(n);
        for (int i = 2; i < n; ++i)
            F[i] = F[i - 2] + F[i - 1];
        solution A(a), B(b), C, De;
        C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
        De.x = A.x + B.x - C.x;

        C.fit_fun(ff, ud1, ud2);
        De.fit_fun(ff, ud1, ud2);
        for (int i = 0; i <= n - 3; ++i)
        {
            if (C.y < De.y)
                B = De;
            else
                A = C;
            C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
            De.x = A.x + B.x - C.x;
            C.fit_fun(ff, ud1, ud2);
            De.fit_fun(ff, ud1, ud2);
            Xopt.ud.add_row((B.x - A.x)());
        }
        Xopt = C;
        Xopt.flag = 0;
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution fib(...):\n" + ex_info);
    }
}

solution
lag(matrix(*ff)(matrix, matrix, matrix), const double a, const double b, double epsilon, double gamma, int Nmax, matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;
        Xopt.ud = b - a;
        solution A(a), B(b), C, De, D_old(a);
        C.x = (a + b) / 2;
        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);
        double l, m;
        while (true) {
            l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y * (pow(C.x) - pow(A.x)) + C.y * (pow(A.x) - pow(B.x)));
            m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));
            if (m <= 0) {
                Xopt = D_old;
                Xopt.flag = 2;
                return Xopt;
            }
            De.x = 0.5 * l / m;
            De.fit_fun(ff, ud1, ud2);
            if (A.x <= De.x && De.x <= C.x) {
                if (De.y < C.y) {
                    B = C;
                    C = De;
                }
                else {
                    A = De;
                }
            }
            else if (C.x <= De.x && De.x <= B.x) {
                if (De.y < C.y) {
                    A = C;
                    C = De;
                }
                else {
                    B = De;
                }
            }
            else {
                Xopt = D_old;
                Xopt.flag = 2;
                return Xopt;
            }
            Xopt.ud.add_row((B.x - A.x)());
            if (B.x - A.x < epsilon || abs(De.x() - D_old.x()) < gamma) {
                Xopt = De;
                Xopt.flag = 0;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt = De;
                Xopt.flag = 1;
                break;
            }
            D_old = De;
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution lag(...):\n" + ex_info);
    }
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {

    const auto& probuj = [&](solution X, double s) {
        int dim = get_len(X.x);
        std::vector<matrix> e;
        std::vector<std::vector<double>> wersory(dim);
        for (int i = 0; i < wersory.size(); ++i) {
            wersory[i].resize(dim, 0);
            wersory[i][i] = 1.;
            e.emplace_back(dim, wersory[i].data());
        }

        for (int j = 0; j < get_len(X.x); ++j) {
            solution X1(X.x - s * e[j]), X2 = X, X3(X.x + s * e[j]);
            X1.fit_fun(ff, ud1, ud2);
            X3.fit_fun(ff, ud1, ud2);
            if (X3.y < X2.y) {
                X = X3;
            }
            else if (X1.y < X2.y) {
                X = X1;
            }
        }
        return X;
        };

    solution X(x0);
    solution Xb;
    X.fit_fun(ff, ud1, ud2);
    X.ud = trans(matrix(2, std::vector<double>{X.x(0), X.x(1)}.data()));

    while (true) {
        Xb = X;
        X = probuj(Xb, s);
        if (X.y < Xb.y) {
            while (true) {
                solution Xb_tmp = Xb;
                Xb = X;
                X.x = 2 * Xb.x - Xb_tmp.x;
                X = probuj(Xb, s);
                if (solution::f_calls > Nmax) {
                    throw std::runtime_error("Fcalls greater than Nmax, fcalls = "
                        + std::to_string(solution::f_calls));
                }
                if (!(X.y < Xb.y)) break;
            }
            X = Xb;
        }
        else {
            s *= alpha;
        }
        if (solution::f_calls > Nmax) {
            throw std::runtime_error("Fcalls greater than Nmax, fcalls = "
                + std::to_string(solution::f_calls));
        }
        X.ud.add_row(trans(matrix(2, std::vector<double>{X.x(0), X.x(1)}.data())));
        if (s < epsilon) break;
    }
    return Xb;
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    int i = 0;
    int dim = get_len(x0);
    matrix d = ident_mat(dim);
    matrix lambda(dim, std::vector<double>(dim, 0).data());
    matrix p(dim, std::vector<double>(dim, 0).data());

    solution X(x0);
    X.fit_fun(ff, ud1, ud2);
    solution Xb = X;
    ud1 = trans(matrix(2, std::vector<double>{X.x(0), X.x(1)}.data()));
    while (true) {
        s0.add_col();
        lambda.add_col();
        p.add_col();
        for (int j = 0; j < dim; ++j) {
            solution Xb_tmp = Xb;
            Xb_tmp.x = Xb_tmp.x + s0(j, i) * d[j];
            Xb_tmp.fit_fun(ff, ud1, ud2);
            if (Xb_tmp.y < Xb.y) {
                Xb = Xb_tmp;
                lambda(j, i + 1) = lambda(j, i) + s0(j, i);
                s0(j, i + 1) = alpha * s0(j, i);
                p(j, i + 1) = p(j, i); // bierzemy z poprzedniej iteracji?
            }
            else
            {
                s0(j, i + 1) = -beta * s0(j, i);
                p(j, i + 1) = p(j, i) + 1;
                lambda(j, i + 1) = lambda(j, i); // bierzemy z poprzedniej iteracji?
            }
        }
        ++i;
        X = Xb;
        bool czy_zmienic_baze = true;
        for (int j = 0; j < dim; ++j) {
            if (!(lambda(j, i) != 0 && p(j, i) != 0)) {
                czy_zmienic_baze = false;
                break;
            }
        }

        if (czy_zmienic_baze) {
            matrix tmp(dim, dim, 0.);
            for (int j = 0; j < dim; ++j) {
                for (int k = 0; k <= j; ++k) {
                    tmp(j, k) = lambda(j, i);
                }
            }
            matrix q = d * tmp;
            for (int j = 0; j < dim; ++j) {
                matrix sum;
                for (int k = 0; k < j; ++k) {
                    sum = sum + (trans(get_col(q, j)) * d[k]) * d[k];
                }
                matrix v = get_col(q, j) - sum;
                v = v / norm(v);
                if (!j) {
                    d = matrix(v);
                }
                else {
                    d.add_col(v);
                }
            }
            lambda = get_col(lambda, 0);
            p = get_col(p, 0);
            s0 = get_col(s0, 0);
            i = 0;
        }

        ud1.add_row(trans(matrix(2, std::vector<double>{X.x(0), X.x(1)}.data())));
        if (solution::f_calls > Nmax) {
            throw std::runtime_error("Fcalls greater than Nmax, fcalls = "
                + std::to_string(solution::f_calls));
        }
        double max{};
        for (int j = 0; j < dim; ++j) {
            if (fabs(s0(j, i)) > max) {
                max = fabs(s0(j, i));
            }
        }
        if (max < epsilon) break;
    }
    return X;
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{

    solution Xopt;
    int n = get_len(x0);
    matrix De = ident_mat(n);
    int N = n + 1;
    solution* S = new solution[N];
    S[0].x = x0;
    S[0].fit_fun(ff, ud1, ud2);
    for (int i = 1; i < N; ++i)
    {
        S[i].x = S[0].x + s * De[i - 1];
        S[i].fit_fun(ff, ud1, ud2);
    }
    solution PR, PE, PN;
    matrix pc;
    int i_min, i_max;
    while (true)
    {
        i_min = i_max = 0;
        for (int i = 1; i < N; ++i)
        {
            if (S[i_min].y > S[i].y) {
                i_min = i;
            }
            if (S[i_max].y < S[i].y) {
                i_max = i;
            }

        }
        pc = matrix(n, 1);
        for (int i = 0; i < N; ++i)
            if (i != i_max) {
                pc = pc + S[i].x;
            }
        pc = pc / (N - 1);
        PR.x = pc + alpha * (pc - S[i_max].x);
        PR.fit_fun(ff, ud1, ud2);
        if (PR.y < S[i_max].y && S[i_min].y <= PR.y) {
            S[i_max] = PR;
        }
        else if (PR.y < S[i_min].y)
        {
            PE.x = pc + gamma * (PR.x - pc);
            PE.fit_fun(ff, ud1, ud2);
            if (PE.y < PR.y) {
                S[i_max] = PE;
            }
            else {
                S[i_max] = PR;
            }
        }
        else
        {
            PN.x = pc + beta * (S[i_max].x - pc);
            PN.fit_fun(ff, ud1, ud2);
            if (PN.y < S[i_max].y) {
                S[i_max] = PN;
            }
            else
            {
                for (int i = 0; i < N; ++i) {
                    if (i != i_min) {
                        S[i].x = delta * (S[i].x + S[i_min].x);
                        S[i].fit_fun(ff, ud1, ud2);
                    }
                }
            }
        }
        double max_s = norm(S[0].x - S[i_min].x);
        for (int i = 1; i < N; ++i) {
            if (max_s < norm(S[i].x - S[i_min].x)) {
                max_s = norm(S[i].x - S[i_min].x);
            }
        }
        if (max_s < epsilon || solution::f_calls > Nmax) {
            Xopt = S[i_min];
            break;
        }
    }
    return Xopt;
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c0, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
    solution X(x0), X1;
    matrix c(2, new double[2] { c0, dc });
    while (true)
    {
        X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
        if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax) {
            return X1;
        }
        c(0) = dc * c(0);
        X = X1;
    }
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = trans(x0);
        int n = get_len(x0);
        solution X0, X1;
        X0.x = x0;
        matrix d(n, 1), P(n, 2);
        solution h;
        double* ab{};

        while (true)
        {
            X0.grad(gf, ud1, ud2);
            d = -X0.g;
            if (h0 < 0)
            {
                P.set_col(X0.x, 0);
                P.set_col(d, 1);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
                X1.x = X0.x + h.x * d;
            }
            else
                X1.x = X0.x + h0 * d;

            Xopt.ud.add_row(trans(X1.x));
            if (norm(X0.x - X1.x) < epsilon)
            {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if (std::max(solution::f_calls, solution::g_calls) > Nmax)
            {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }
            X0 = X1;
            //            std::cout << "X0.x = " << X0.x << std::endl;
        }
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = trans(x0);
        int n = get_len(x0);
        solution X0, X1;
        X0.x = x0;
        matrix d(n, 1), P(n, 2);
        solution h;
        double* ab{};
        double beta;
        d = -X0.grad(gf);
        while (true)
        {
            if (h0 < 0)
            {
                P.set_col(X0.x, 0);
                P.set_col(d, 1);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
                X1.x = X0.x + h.x * d;
            }
            else
                X1.x = X0.x + h0 * d;

            Xopt.ud.add_row(trans(X1.x));

            if (norm(X1.x - X0.x) < epsilon)
            {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if (std::max(solution::f_calls, solution::g_calls) > Nmax)
            {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }
            X1.grad(gf);
            beta = pow(norm(X1.g), 2) / pow(norm(X0.g), 2);
            d = -X1.g + beta * d;
            X0 = X1;
            //            std::cout << "X0.x = " << X0.x << std::endl;
        }
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
    matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = trans(x0);
        int n = get_len(x0);
        solution X0, X1;
        X0.x = x0;
        matrix d(n, 1), P(n, 2);
        solution h;
        double* ab{};
        while (true)
        {
            X0.grad(gf);
            X0.hess(Hf);
            d = -inv(X0.H) * X0.g;
            if (h0 < 0)
            {
                P.set_col(X0.x, 0);
                P.set_col(d, 1);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
                X1.x = X0.x + h.x * d;
            }
            else
                X1.x = X0.x + h0 * d;

            Xopt.ud.add_row(trans(X1.x));

            if (norm(X0.x - X1.x) < epsilon)
            {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if (std::max(solution::f_calls, solution::g_calls) > Nmax)
            {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }
            //            std::cout << "X0.x = " << X0.x << std::endl;
            X0 = X1;
        }
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        double alfa = (sqrt(5) - 1) / 2;
        solution A, B, C, De;
        A.x = a;
        B.x = b;
        C.x = B.x - alfa * (B.x - A.x);
        C.fit_fun(ff, ud1, ud2);
        De.x = A.x + alfa * (B.x - A.x);
        De.fit_fun(ff, ud1, ud2);
        while (true)
        {
            if (C.y < De.y)
            {
                B = De;
                De = C;
                C.x = B.x - alfa * (B.x - A.x);
                C.fit_fun(ff, ud1, ud2);
            }
            else
            {
                A = C;
                C = De;
                De.x = A.x + alfa * (B.x - A.x);
                De.fit_fun(ff, ud1, ud2);
            }
            if (B.x - A.x < epsilon)
            {
                Xopt.x = (A.x + B.x) / 2;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if (std::max(solution::f_calls, solution::g_calls) > Nmax)
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
        matrix De = ident_mat(n), A(n, 2);
        solution X, P, h;
        X.x = x0;
        double* ab;
        while (true)
        {
            P = X;
            for (int i = 0; i < n; ++i)
            {
                A.set_col(P.x, 0);
                A.set_col(De[i], 1);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
                P.x = P.x + h.x * De[i];
            }
            if (norm(X.x) < epsilon)
            {
                Xopt = P;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if (solution::f_calls > Nmax)
            {
                Xopt = P;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }
            for (int i = 0; i < n - 1; ++i) De.set_col(De[i + 1], i);
            De.set_col(P.x - X.x, n - 1);
            A.set_col(P.x, 0);
            A.set_col(De[n - 1], 1);
            ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
            h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
            /// ?????
            X.x = P.x;
        }
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        solution* P = new solution[mi + lambda];
        solution* Pm = new solution[mi];
        matrix IFF(mi, 1), temp(N, 2);
        double r, s, s_IFF;
        double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
        int j_min;
        for (int i = 0; i < mi; ++i)
        {
            P[i].x = matrix(N, 2);
            for (int j = 0; j < N; ++j)
            {
                P[i].x(j, 0) = (limits(j, 1) - limits(j, 0)) * m2d(rand_mat()) + limits(j, 0);
                P[i].x(j, 1) = sigma0(j);
            }
            if (P[i].fit_fun(ff, ud1, ud2) < epsilon)
            {
                Xopt = P[i];
                Xopt.flag = 0;
                delete[]P;
                delete[]Pm;
                return Xopt;
            }
        }
        while (true)
        {
            s_IFF = 0;
            for (int i = 0; i < mi; ++i)
            {
                IFF(i) = 1 / m2d(P[i].y);
                s_IFF += IFF(i);
            }
            for (int i = 0; i < lambda; ++i)
            {
                r = s_IFF * m2d(rand_mat());
                s = 0;
                for (int j = 0; j < mi; ++j)
                {
                    s += IFF(j);
                    if (r <= s)
                    {
                        P[mi + i] = P[j];
                        break;
                    }
                }
            }
            for (int i = 0; i < lambda; ++i)
            {
                r = m2d(randn_mat());
                for (int j = 0; j < N; ++j)
                {
                    P[mi + i].x(j, 1) *= exp(tau1 * r + tau * m2d(randn_mat()));
                    P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * m2d(randn_mat());
                }
            }
            for (int i = 0; i < lambda; i += 2)
            {
                r = m2d(rand_mat());
                temp = P[mi + i].x;
                P[mi + i].x = r * P[mi + i].x + (1 - r) * P[mi + i + 1].x;
                P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r) * temp;
            }
            for (int i = 0; i < lambda; ++i)
            {
                if (P[mi + i].fit_fun(ff, ud1, ud2) < epsilon)
                {
                    Xopt = P[mi + i];
                    Xopt.flag = 0;
                    delete[]P;
                    delete[]Pm;
                    return Xopt;
                }
            }
            for (int i = 0; i < mi; ++i)
            {
                j_min = 0;
                for (int j = 1; j < mi + lambda; ++j)
                    if (P[j_min].y > P[j].y)
                        j_min = j;
                Pm[i] = P[j_min];
                P[j_min].y = 1e10;
            }
            for (int i = 0; i < mi; ++i)
                P[i] = Pm[i];
            if (solution::f_calls > Nmax)
            {
                Xopt = P[0];
                Xopt.flag = 1;
                break;
            }
        }

        delete[]P;
        delete[]Pm;
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution EA(...):\n" + ex_info);
    }
}