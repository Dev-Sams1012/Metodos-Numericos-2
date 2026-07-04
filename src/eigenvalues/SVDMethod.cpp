#include "eigenvalues/SVDMethod.hpp"
#include "eigenvalues/JacobiMethod.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace nm
{

    SVDResult SVDSolver::solve(const Matrix &A) const
    {
        int m = A.rows();
        int n = A.cols();

        Matrix At = A.transpose();
        Matrix B = At * A;

        JacobiMethod jacobi;
        JacobiResult jr = jacobi.solve(B, tol_ * tol_, 50000);

        int nEig = (int)jr.eigenvalues.size();

        std::vector<int> idx(nEig);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b)
                  { return jr.eigenvalues[a] > jr.eigenvalues[b]; });

        std::vector<double> sigmas;
        sigmas.reserve(nEig);
        int rank = 0;
        for (int k : idx)
        {
            double lam = jr.eigenvalues[k];
            double sig = (lam > tol_ * tol_) ? std::sqrt(std::max(0.0, lam)) : 0.0;
            sigmas.push_back(sig);
            if (sig > tol_)
                ++rank;
        }

        Matrix V(n, rank, 0.0);
        for (int col = 0; col < rank; ++col)
        {
            int srcCol = idx[col];
            for (int row = 0; row < n; ++row)
                V(row, col) = jr.eigenvectors(row, srcCol);
        }

        Matrix U(m, rank, 0.0);
        for (int col = 0; col < rank; ++col)
        {
            double sig = sigmas[col];

            std::vector<double> vi(n);
            for (int row = 0; row < n; ++row)
                vi[row] = V(row, col);

            std::vector<double> ui = A.multiply(vi);
            for (int row = 0; row < m; ++row)
                U(row, col) = ui[row] / sig;
        }

        std::vector<double> svs(sigmas.begin(), sigmas.begin() + rank);

        return {U, svs, V, rank};
    }

}