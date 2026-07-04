#pragma once

#include "ivp/IVPSolver.hpp"

namespace nm
{

    class RungeKutta2 : public IVPSolver
    {
    public:
        RungeKutta2() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override { return "Runge-Kutta 2ª ordem  O(Δt²)  (Heun)"; }

        std::vector<double> step(const SystemParser &F,
                                 const std::vector<double> &y,
                                 double t, double dt) const;
    };

    class RungeKutta3 : public IVPSolver
    {
    public:
        RungeKutta3() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override { return "Runge-Kutta 3ª ordem  O(Δt³)  (Kutta/Simpson)"; }

        std::vector<double> step(const SystemParser &F,
                                 const std::vector<double> &y,
                                 double t, double dt) const;
    };

    class RungeKutta4 : public IVPSolver
    {
    public:
        RungeKutta4() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override { return "Runge-Kutta 4ª ordem  O(Δt⁴)  (Clássico)"; }

        std::vector<double> step(const SystemParser &F,
                                 const std::vector<double> &y,
                                 double t, double dt) const;
    };

}