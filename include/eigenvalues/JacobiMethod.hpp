#pragma once

#include "eigenvalues/EigenSolver.hpp"
#include <vector>

namespace nm
{

    struct JacobiResult
    {
        std::vector<double> eigenvalues;
        Matrix eigenvectors;
        int iterations;
        bool converged;
    };

    class JacobiMethod
    {
    public:
        JacobiMethod() = default;

        JacobiResult solve(const Matrix &A, double tol = 1e-10,
                           int maxIter = 10000) const;

        std::string name() const { return "Método de Jacobi"; }
    };

}