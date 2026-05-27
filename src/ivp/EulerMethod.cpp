// =============================================================================
// EulerMethod.cpp
// =============================================================================
#include "ivp/EulerMethod.hpp"
#include <cmath>

namespace nm
{

    // ============================================================
    //  Euler Explícito
    // ============================================================
    IVPResult EulerExplicit::solve(const SystemParser &F,
                                   const std::vector<double> &y0,
                                   double t0, double tEnd, double dt,
                                   int maxSteps) const
    {
        IVPResult res;
        res.completed = false;

        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        int step = 0;
        while (t < tEnd - 1e-12 * std::abs(dt) && step < maxSteps)
        {
            double dtEff = std::min(dt, tEnd - t);

            // S_{i+1} = S_i + dt * F(S_i, t_i)
            auto f = F.eval(y, t);
            y = axpy(y, dtEff, f);
            t += dtEff;
            ++step;

            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = step;
        res.completed = (step < maxSteps);
        return res;
    }

    // ============================================================
    //  Euler Implícito  (iteração de ponto fixo)
    // ============================================================
    IVPResult EulerImplicit::solve(const SystemParser &F,
                                   const std::vector<double> &y0,
                                   double t0, double tEnd, double dt,
                                   int maxSteps) const
    {
        IVPResult res;
        res.completed = false;

        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        int step = 0;
        while (t < tEnd - 1e-12 * std::abs(dt) && step < maxSteps)
        {
            double dtEff = std::min(dt, tEnd - t);
            double tNext = t + dtEff;

            // Chute inicial: Euler explícito
            auto fi = F.eval(y, t);
            auto yNew = axpy(y, dtEff, fi);

            // Iteração de ponto fixo
            for (int iter = 0; iter < maxIter_; ++iter)
            {
                auto yPrev = yNew;
                auto fNew = F.eval(yNew, tNext);
                yNew = axpy(y, dtEff, fNew);

                if (relDiff(yNew, yPrev) < tol_)
                    break;
            }

            y = yNew;
            t = tNext;
            ++step;

            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = step;
        res.completed = (step < maxSteps);
        return res;
    }

} // namespace nm