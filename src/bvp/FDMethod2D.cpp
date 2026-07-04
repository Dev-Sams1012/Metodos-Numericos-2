#include "bvp/FDMethod2D.hpp"
#include <stdexcept>
#include <cmath>

namespace nm
{

    BVPResult FDMethod2D::solve(std::function<double(double, double)> f,
                                double xa, double xb,
                                double ya, double yb,
                                int nx, int ny,
                                double fVal) const
    {
        if (ny <= 0)
            ny = nx;
        if (nx < 2 || ny < 2)
            throw std::invalid_argument("FDMethod2D: nx e ny devem ser >= 2.");
        if (xb <= xa || yb <= ya)
            throw std::invalid_argument("FDMethod2D: domínio inválido.");

        const double hx = (xb - xa) / static_cast<double>(nx);
        const double hy = (yb - ya) / static_cast<double>(ny);
        const double hx2 = hx * hx;
        const double hy2 = hy * hy;

        const int mx = nx - 1;
        const int my = ny - 1;
        const int M = mx * my;

        const double cW = 1.0 / hx2;
        const double cE = 1.0 / hx2;
        const double cS = 1.0 / hy2;
        const double cN = 1.0 / hy2;
        const double cC = -2.0 / hx2 - 2.0 / hy2;

        auto K = [&](int i, int j) -> int
        {
            return (j - 1) * mx + (i - 1);
        };

        BVPResult result;
        result.n = nx;
        result.A = Matrix(M, M, 0.0);
        result.rhs.assign(M, 0.0);

        for (int j = 1; j <= my; ++j)
        {
            double yj = ya + j * hy;
            for (int i = 1; i <= mx; ++i)
            {
                double xi = xa + i * hx;
                int k = K(i, j);

                result.A(k, k) = cC;

                if (i > 1)
                    result.A(k, K(i - 1, j)) += cW;

                if (i < mx)
                    result.A(k, K(i + 1, j)) += cE;

                if (j > 1)
                    result.A(k, K(i, j - 1)) += cS;

                if (j < my)
                    result.A(k, K(i, j + 1)) += cN;

                double fij = f ? f(xi, yj) : fVal;
                result.rhs[k] = fij;
            }
        }

        result.x.resize(M);
        result.y.resize(M);
        for (int j = 1; j <= my; ++j)
            for (int i = 1; i <= mx; ++i)
            {
                int k = K(i, j);
                result.x[k] = xa + i * hx;
                result.y[k] = ya + j * hy;
            }

        try
        {
            result.u = Matrix::solve(result.A, result.rhs);
            result.solved = true;
        }
        catch (const std::runtime_error &)
        {
            result.solved = false;
        }

        return result;
    }

}