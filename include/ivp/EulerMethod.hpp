#pragma once

#include "ivp/IVPSolver.hpp"

namespace nm
{

    class EulerExplicit : public IVPSolver
    {
    public:
        EulerExplicit() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override { return "Euler Explícito  O(Δt)"; }
    };

    class EulerImplicit : public IVPSolver
    {
    public:
        explicit EulerImplicit(double tol = 1e-10, int maxIter = 100)
            : tol_(tol), maxIter_(maxIter) {}

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override { return "Euler Implícito  O(Δt)  (ponto fixo)"; }

    private:
        double tol_;
        int maxIter_;
    };

}