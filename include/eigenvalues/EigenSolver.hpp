#pragma once

#include "Matrix.hpp"
#include <vector>
#include <string>
#include <stdexcept>

namespace nm
{

    struct EigenResult
    {
        double eigenvalue;
        std::vector<double> eigenvector;
        int iterations;
        bool converged;
    };

    class EigenSolver
    {
    public:
        virtual ~EigenSolver() = default;

        EigenSolver(const EigenSolver &) = delete;
        EigenSolver &operator=(const EigenSolver &) = delete;

        virtual EigenResult solve(const Matrix &A, double tol = 1e-10,
                                  int maxIter = 1000) const = 0;

        virtual std::string name() const = 0;

    protected:
        EigenSolver() = default;
    };

}