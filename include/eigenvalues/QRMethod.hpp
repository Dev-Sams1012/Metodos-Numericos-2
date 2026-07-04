#pragma once

#include "eigenvalues/EigenSolver.hpp"
#include <vector>

namespace nm
{

    struct QRResult
    {
        std::vector<double> eigenvalues;
        Matrix eigenvectors;
        int iterations;
        bool converged;
    };

    class QRMethod
    {
    public:
        explicit QRMethod(bool useHouseholder = true)
            : useHouseholder_(useHouseholder) {}

        QRResult solve(const Matrix &A, double tol = 1e-10,
                       int maxIter = 10000) const;

        std::string name() const { return "Método QR (decomposição iterativa)"; }

    private:
        bool useHouseholder_;

        void qrDecompose(const Matrix &A, Matrix &Q, Matrix &R) const;

        bool isConverged(const Matrix &A, double tol) const;
    };

}