#include "ivp/RungeKutta.hpp"
#include <cmath>

namespace nm
{

    static std::vector<double> add3(const std::vector<double> &a,
                                    double c1, const std::vector<double> &b1,
                                    double c2, const std::vector<double> &b2)
    {
        int n = (int)a.size();
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i)
            r[i] = a[i] + c1 * b1[i] + c2 * b2[i];
        return r;
    }

    static std::vector<double> add5(const std::vector<double> &a,
                                    double c1, const std::vector<double> &k1,
                                    double c2, const std::vector<double> &k2,
                                    double c3, const std::vector<double> &k3,
                                    double c4, const std::vector<double> &k4)
    {
        int n = (int)a.size();
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i)
            r[i] = a[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i] + c4 * k4[i];
        return r;
    }

    std::vector<double> RungeKutta2::step(const SystemParser &F,
                                          const std::vector<double> &y,
                                          double t, double dt) const
    {
        auto k1 = F.eval(y, t);
        auto yEst = IVPSolver::axpy(y, dt, k1);
        auto k2 = F.eval(yEst, t + dt);
        return add3(y, dt / 2.0, k1, dt / 2.0, k2);
    }

    IVPResult RungeKutta2::solve(const SystemParser &F,
                                 const std::vector<double> &y0,
                                 double t0, double tEnd, double dt,
                                 int maxSteps) const
    {
        IVPResult res;
        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        int s = 0;
        while (t < tEnd - 1e-12 * std::abs(dt) && s < maxSteps)
        {
            double h = std::min(dt, tEnd - t);
            y = step(F, y, t, h);
            t += h;
            ++s;
            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = s;
        res.completed = (s < maxSteps);
        return res;
    }

    std::vector<double> RungeKutta3::step(const SystemParser &F,
                                          const std::vector<double> &y,
                                          double t, double dt) const
    {

        auto k1 = F.eval(y, t);

        auto yHalf = IVPSolver::axpy(y, dt / 2.0, k1);
        auto k2 = F.eval(yHalf, t + dt / 2.0);

        auto yFull = add3(y, -dt, k1, 2.0 * dt, k2);
        auto k3 = F.eval(yFull, t + dt);

        int n = (int)y.size();
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i)
            r[i] = y[i] + (dt / 6.0) * (k1[i] + 4.0 * k2[i] + k3[i]);
        return r;
    }

    IVPResult RungeKutta3::solve(const SystemParser &F,
                                 const std::vector<double> &y0,
                                 double t0, double tEnd, double dt,
                                 int maxSteps) const
    {
        IVPResult res;
        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        int s = 0;
        while (t < tEnd - 1e-12 * std::abs(dt) && s < maxSteps)
        {
            double h = std::min(dt, tEnd - t);
            y = step(F, y, t, h);
            t += h;
            ++s;
            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = s;
        res.completed = (s < maxSteps);
        return res;
    }

    std::vector<double> RungeKutta4::step(const SystemParser &F,
                                          const std::vector<double> &y,
                                          double t, double dt) const
    {
        auto k1 = F.eval(y, t);
        auto k2 = F.eval(IVPSolver::axpy(y, dt / 2.0, k1), t + dt / 2.0);
        auto k3 = F.eval(IVPSolver::axpy(y, dt / 2.0, k2), t + dt / 2.0);
        auto k4 = F.eval(IVPSolver::axpy(y, dt, k3), t + dt);
        return add5(y, dt / 6.0, k1, dt / 3.0, k2, dt / 3.0, k3, dt / 6.0, k4);
    }

    IVPResult RungeKutta4::solve(const SystemParser &F,
                                 const std::vector<double> &y0,
                                 double t0, double tEnd, double dt,
                                 int maxSteps) const
    {
        IVPResult res;
        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        int s = 0;
        while (t < tEnd - 1e-12 * std::abs(dt) && s < maxSteps)
        {
            double h = std::min(dt, tEnd - t);
            y = step(F, y, t, h);
            t += h;
            ++s;
            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = s;
        res.completed = (s < maxSteps);
        return res;
    }

}