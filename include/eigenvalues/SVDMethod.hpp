#pragma once

#include "Matrix.hpp"
#include <vector>
#include <string>

namespace nm
{

    struct SVDResult
    {
        Matrix U;
        std::vector<double> singularValues;
        Matrix V;
        int rank;
    };

    class SVDSolver
    {
    public:
        explicit SVDSolver(double tol = 1e-10) : tol_(tol) {}

        SVDResult solve(const Matrix &A) const;

        std::string name() const
        {
            return "SVD — Decomposição em Valores Singulares (via autovalores de AᵀA)";
        }

    private:
        double tol_;
    };

}