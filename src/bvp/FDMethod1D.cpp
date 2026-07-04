#include "bvp/FDMethod1D.hpp"
#include <stdexcept>
#include <cmath>

namespace nm
{

    BVPResult FDMethod1D::solve(const std::shared_ptr<FunctionParser> &f,
                                double a, double b,
                                const BoundaryCondition &bcA,
                                const BoundaryCondition &bcB,
                                int n) const
    {
        if (n < 2)
            throw std::invalid_argument("FDMethod1D: n deve ser >= 2.");
        if (b <= a)
            throw std::invalid_argument("FDMethod1D: b deve ser > a.");
        if (std::abs(alpha_) < 1e-14)
            throw std::invalid_argument("FDMethod1D: alpha não pode ser zero.");

        const double h = (b - a) / static_cast<double>(n);
        const double h2 = h * h;

        BVPResult result;
        result.n = n;
        result.x.resize(n + 1);
        result.u.resize(n + 1, 0.0);
        for (int i = 0; i <= n; ++i)
            result.x[i] = a + i * h;

        const double cL = alpha_ / h2 - beta_ / (2.0 * h);
        const double cC = -2.0 * alpha_ / h2 + gamma_;
        const double cR = alpha_ / h2 + beta_ / (2.0 * h);

        const int iFirst = (bcA.type == BCType::DIRICHLET) ? 1 : 0;
        const int iLast = n - 1;
        const int m = iLast - iFirst + 1;

        result.A = Matrix(m, m, 0.0);
        result.rhs.assign(m, 0.0);

        for (int i = iFirst; i <= iLast; ++i)
        {

            const int row = i - iFirst;
            const double xi = result.x[i];
            const double fi = (f && f->isReady()) ? f->eval(xi) : 0.0;

            result.rhs[row] = fi;

            if (i == n - 1 && bcB.type == BCType::NEUMANN)
            {

                if (i - 1 >= iFirst)
                    result.A(row, row - 1) = cL;
                else
                    result.rhs[row] -= cL * bcA.value;

                result.A(row, row) = cC + cR;

                result.rhs[row] -= cR * h * bcB.value;

                continue;
            }

            if (i - 1 < iFirst)
            {
                result.rhs[row] -= cL * bcA.value;
            }
            else
            {
                result.A(row, row - 1) = cL;
            }

            result.A(row, row) = cC;

            if (i + 1 > iLast)
            {

                if (bcB.type == BCType::DIRICHLET)
                    result.rhs[row] -= cR * bcB.value;
            }
            else
            {
                result.A(row, row + 1) = cR;
            }
        }

        try
        {
            std::vector<double> sol = Matrix::solve(result.A, result.rhs);

            for (int i = iFirst; i <= iLast; ++i)
                result.u[i] = sol[i - iFirst];

            if (bcA.type == BCType::DIRICHLET)
                result.u[0] = bcA.value;

            if (bcB.type == BCType::DIRICHLET)
                result.u[n] = bcB.value;

            if (bcB.type == BCType::NEUMANN)
                result.u[n] = result.u[n - 1] + h * bcB.value;

            result.solved = true;
        }
        catch (const std::runtime_error &)
        {
            result.solved = false;
        }

        return result;
    }

}