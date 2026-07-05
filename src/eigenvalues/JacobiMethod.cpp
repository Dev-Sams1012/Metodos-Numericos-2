#include "eigenvalues/JacobiMethod.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace nm
{

    JacobiResult JacobiMethod::solve(const Matrix &A, double tol, int maxIter) const
    {
        if (!A.isSquare())
            throw std::invalid_argument("Jacobi: matriz deve ser quadrada.");

        int n = A.rows();

        Matrix S = A;
        Matrix E = Matrix::identity(n);

        std::vector<double> e(n);
        std::vector<int> ind(n);
        std::vector<bool> changed(n, true);

        int state = n;
        int iter = 0;
        bool converged = false;

        auto maxind = [&](int k) -> int
        {
            if (k >= n - 1)
                return k;
            int m = k + 1;
            for (int i = k + 2; i < n; ++i)
            {
                if (std::abs(S(k, i)) > std::abs(S(k, m)))
                {
                    m = i;
                }
            }
            return m;
        };

        auto update = [&](int k, double t)
        {
            double y = e[k];
            e[k] = y + t;
            if (changed[k] && (y == e[k]))
            {
                changed[k] = false;
                state--;
            }
            else if (!changed[k] && (y != e[k]))
            {
                changed[k] = true;
                state++;
            }
        };

        auto apply_rotate = [&](double &skl, double &sij, double c, double s)
        {
            double temp_kl = c * skl - s * sij;
            double temp_ij = s * skl + c * sij;
            skl = temp_kl;
            sij = temp_ij;
        };

        for (int k = 0; k < n; ++k)
        {
            ind[k] = maxind(k);
            e[k] = S(k, k);
            changed[k] = true;
        }

        while (state != 0 && iter < maxIter)
        {

            int m = 0;
            for (int k = 1; k < n - 1; ++k)
            {
                if (std::abs(S(k, ind[k])) > std::abs(S(m, ind[m])))
                {
                    m = k;
                }
            }

            int k = m;
            int l = ind[m];
            double p = S(k, l);

            if (std::abs(p) <= tol)
            {
                converged = true;
                break;
            }

            double y = (e[l] - e[k]) / 2.0;
            double d = std::abs(y) + std::sqrt(p * p + y * y);

            if (d == 0.0)
                break;

            double r = std::sqrt(p * p + d * d);
            double c = d / r;
            double s = p / r;
            double t = (p * p) / d;

            if (y < 0.0)
            {
                s = -s;
                t = -t;
            }

            S(k, l) = 0.0;
            update(k, -t);
            update(l, t);

            for (int i = 0; i < k; ++i)
                apply_rotate(S(i, k), S(i, l), c, s);
            for (int i = k + 1; i < l; ++i)
                apply_rotate(S(k, i), S(i, l), c, s);
            for (int i = l + 1; i < n; ++i)
                apply_rotate(S(k, i), S(l, i), c, s);

            for (int i = 0; i < n; ++i)
            {
                apply_rotate(E(i, k), E(i, l), c, s);
            }

            for (int i = 0; i < n - 1; ++i)
            {
                ind[i] = maxind(i);
            }

            iter++;
        }

        if (state == 0)
        {
            converged = true;
        }

        return {e, E, iter, converged};
    }

}