#pragma once

#include "ivp/IVPSolver.hpp"

namespace nm
{

    class AdamsBashforth2 : public IVPSolver
    {
    public:
        AdamsBashforth2() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override
        {
            return "Adams-Bashforth-Moulton 2ª ordem  O(Δt²)  (PC, init: RK2)";
        }
    };

    class AdamsBashforth3 : public IVPSolver
    {
    public:
        AdamsBashforth3() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override
        {
            return "Adams-Bashforth-Moulton 3ª ordem  O(Δt³)  (PC, init: RK3)";
        }
    };

    class AdamsBashforth4 : public IVPSolver
    {
    public:
        AdamsBashforth4() = default;

        IVPResult solve(const SystemParser &F,
                        const std::vector<double> &y0,
                        double t0, double tEnd, double dt,
                        int maxSteps = 200000) const override;

        std::string name() const override
        {
            return "Adams-Bashforth-Moulton 4ª ordem  O(Δt⁴)  (PC, init: RK4)";
        }
    };

}