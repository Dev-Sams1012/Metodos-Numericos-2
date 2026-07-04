#include "eigenvalues/HouseholderMethod.hpp"
#include <cmath>
#include <stdexcept>

namespace nm
{

    HouseholderResult HouseholderMethod::tridiagonalize(const Matrix &A) const
    {
        if (!A.isSquare())
            throw std::invalid_argument("Householder: matriz deve ser quadrada.");
        int n = A.rows();

        Matrix T = A;
        Matrix Q = Matrix::identity(n);

        for (int k = 0; k < n - 2; ++k)
        {

            int len = n - k - 1;
            std::vector<double> x(len);
            for (int i = 0; i < len; ++i)
                x[i] = T(k + 1 + i, k);

            double normX = Matrix::vecNorm(x);
            if (normX < 1e-14)
                continue;

            double alpha = (x[0] >= 0.0) ? -normX : normX;

            std::vector<double> w(len);
            w[0] = x[0] - alpha;
            for (int i = 1; i < len; ++i)
                w[i] = x[i];

            double wTw = Matrix::dot(w, w);
            if (wTw < 1e-14)
                continue;

            for (int j = k; j < n; ++j)
            {

                double p = 0.0;
                for (int i = 0; i < len; ++i)
                    p += w[i] * T(k + 1 + i, j);
                p *= 2.0 / wTw;
                for (int i = 0; i < len; ++i)
                    T(k + 1 + i, j) -= p * w[i];
            }

            for (int i = 0; i < n; ++i)
            {

                double q = 0.0;
                for (int j = 0; j < len; ++j)
                    q += T(i, k + 1 + j) * w[j];
                q *= 2.0 / wTw;
                for (int j = 0; j < len; ++j)
                    T(i, k + 1 + j) -= q * w[j];
            }

            for (int i = 0; i < n; ++i)
            {
                double r = 0.0;
                for (int j = 0; j < len; ++j)
                    r += Q(i, k + 1 + j) * w[j];
                r *= 2.0 / wTw;
                for (int j = 0; j < len; ++j)
                    Q(i, k + 1 + j) -= r * w[j];
            }
        }

        return {T, Q, n - 2};
    }

}