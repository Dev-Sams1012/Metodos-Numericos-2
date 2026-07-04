#pragma once

#include "eigenvalues/EigenSolver.hpp"

namespace nm
{

    class PowerRegular : public EigenSolver
    {
    public:
        PowerRegular() = default;
        EigenResult solve(const Matrix &A, double tol = 1e-10,
                          int maxIter = 1000) const override;
        std::string name() const override { return "Método da Potência Regular"; }
    };

    class PowerInverse : public EigenSolver
    {
    public:
        PowerInverse() = default;
        EigenResult solve(const Matrix &A, double tol = 1e-10,
                          int maxIter = 1000) const override;
        std::string name() const override { return "Método da Potência Inverso"; }
    };

    class PowerShifted : public EigenSolver
    {
    public:
        explicit PowerShifted(double mu = 0.0) : mu_(mu) {}

        void setShift(double mu) { mu_ = mu; }
        double shift() const { return mu_; }

        EigenResult solve(const Matrix &A, double tol = 1e-10,
                          int maxIter = 1000) const override;
        std::string name() const override;

    private:
        double mu_;
    };

}