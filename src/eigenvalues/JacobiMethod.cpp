#include "eigenvalues/JacobiMethod.hpp"
#include <cmath>
#include <stdexcept>
#include <limits>

namespace nm
{

    static constexpr double PI = 3.14159265358979323846;

    static double offDiagNorm(const Matrix &A)
    {
        int n = A.rows();
        double s = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    s += A(i, j) * A(i, j);
        return std::sqrt(s);
    }

    static std::pair<int, int> maxOffDiag(const Matrix &A)
    {
        int n = A.rows();
        double maxVal = 0.0;
        int p = 0, q = 1;
        for (int i = 0; i < n - 1; ++i)
            for (int j = i + 1; j < n; ++j)
                if (std::abs(A(i, j)) > maxVal)
                {
                    maxVal = std::abs(A(i, j));
                    p = i;
                    q = j;
                }
        return {p, q};
    }

    JacobiResult JacobiMethod::solve(const Matrix &A, double tol, int maxIter) const
    {
        if (!A.isSquare())
            throw std::invalid_argument("Jacobi: matriz deve ser quadrada.");
        int n = A.rows();

        Matrix S = A;
        Matrix V = Matrix::identity(n);

        int iter = 0;
        bool converged = false;

        for (iter = 1; iter <= maxIter; ++iter)
        {
            if (offDiagNorm(S) < tol)
            {
                converged = true;
                break;
            }

            auto [p, q] = maxOffDiag(S);
            double apq = S(p, q);
            double app = S(p, p);
            double aqq = S(q, q);

            double c, s;
            if (std::abs(app - aqq) < 1e-14)
            {

                c = std::cos(PI / 4.0);
                s = (apq >= 0.0) ? std::sin(PI / 4.0) : -std::sin(PI / 4.0);
            }
            else
            {

                double tau = (aqq - app) / (2.0 * apq);
                double t = (tau >= 0.0)
                               ? 1.0 / (tau + std::sqrt(1.0 + tau * tau))
                               : 1.0 / (tau - std::sqrt(1.0 + tau * tau));
                c = 1.0 / std::sqrt(1.0 + t * t);
                s = t * c;
            }

            for (int r = 0; r < n; ++r)
            {
                if (r == p || r == q)
                    continue;
                double srp = S(r, p);
                double srq = S(r, q);
                S(r, p) = c * srp - s * srq;
                S(p, r) = S(r, p);
                S(r, q) = s * srp + c * srq;
                S(q, r) = S(r, q);
            }

            double new_pp = c * c * app - 2.0 * s * c * apq + s * s * aqq;
            double new_qq = s * s * app + 2.0 * s * c * apq + c * c * aqq;
            S(p, p) = new_pp;
            S(q, q) = new_qq;
            S(p, q) = 0.0;
            S(q, p) = 0.0;

            for (int r = 0; r < n; ++r)
            {
                double vrp = V(r, p);
                double vrq = V(r, q);
                V(r, p) = c * vrp - s * vrq;
                V(r, q) = s * vrp + c * vrq;
            }
        }

        std::vector<double> eigenvalues(n);
        for (int i = 0; i < n; ++i)
            eigenvalues[i] = S(i, i);

        return {eigenvalues, V, iter, converged};
    }

}