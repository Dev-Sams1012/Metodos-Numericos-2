#include "ivp/AdamsBashforth.hpp"
#include "ivp/RungeKutta.hpp"
#include <cmath>

namespace nm
{

    static std::vector<double> linComb(const std::vector<double> &base,
                                       double dt,
                                       std::initializer_list<std::pair<double, const std::vector<double> *>> terms)
    {
        int n = (int)base.size();
        std::vector<double> r = base;
        for (auto &[c, v] : terms)
            for (int i = 0; i < n; ++i)
                r[i] += dt * c * (*v)[i];
        return r;
    }

    IVPResult AdamsBashforth2::solve(const SystemParser &F,
                                     const std::vector<double> &y0,
                                     double t0, double tEnd, double dt,
                                     int maxSteps) const
    {
        IVPResult res;
        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        RungeKutta2 rk2;
        auto fPrev = F.eval(y, t);

        if (t < tEnd - 1e-12)
        {
            double h = std::min(dt, tEnd - t);
            y = rk2.step(F, y, t, h);
            t += h;
            res.t.push_back(t);
            res.y.push_back(y);
        }

        auto fCurr = F.eval(y, t);

        int s = 1;

        while (t < tEnd - 1e-12 * std::abs(dt) && s < maxSteps)
        {
            double h = std::min(dt, tEnd - t);

            auto yPred = linComb(y, h, {{-0.5, &fPrev}, {1.5, &fCurr}});
            auto fPred = F.eval(yPred, t + h);

            auto yNew = linComb(y, h, {{0.5, &fCurr}, {0.5, &fPred}});

            fPrev = fCurr;
            fCurr = F.eval(yNew, t + h);
            y = yNew;
            t += h;
            ++s;

            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = s;
        res.completed = (s < maxSteps);
        return res;
    }

    IVPResult AdamsBashforth3::solve(const SystemParser &F,
                                     const std::vector<double> &y0,
                                     double t0, double tEnd, double dt,
                                     int maxSteps) const
    {
        IVPResult res;
        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        RungeKutta3 rk3;

        auto f0 = F.eval(y, t);

        if (t < tEnd - 1e-12)
        {
            double h = std::min(dt, tEnd - t);
            y = rk3.step(F, y, t, h);
            t += h;
            res.t.push_back(t);
            res.y.push_back(y);
        }
        auto f1 = F.eval(y, t);

        if (t < tEnd - 1e-12)
        {
            double h = std::min(dt, tEnd - t);
            y = rk3.step(F, y, t, h);
            t += h;
            res.t.push_back(t);
            res.y.push_back(y);
        }
        auto f2 = F.eval(y, t);

        int s = 2;

        while (t < tEnd - 1e-12 * std::abs(dt) && s < maxSteps)
        {
            double h = std::min(dt, tEnd - t);

            auto yPred = linComb(y, h / 12.0, {{5.0, &f0}, {-16.0, &f1}, {23.0, &f2}});
            auto fPred = F.eval(yPred, t + h);

            auto yNew = linComb(y, h / 12.0, {{-1.0, &f1}, {8.0, &f2}, {5.0, &fPred}});

            f0 = f1;
            f1 = f2;
            f2 = F.eval(yNew, t + h);
            y = yNew;
            t += h;
            ++s;

            res.t.push_back(t);
            res.y.push_back(y);
        }

        res.steps = s;
        res.completed = (s < maxSteps);
        return res;
    }

    IVPResult AdamsBashforth4::solve(const SystemParser &F,
                                     const std::vector<double> &y0,
                                     double t0, double tEnd, double dt,
                                     int maxSteps) const
    {
        IVPResult res;
        double t = t0;
        std::vector<double> y = y0;

        res.t.push_back(t);
        res.y.push_back(y);

        RungeKutta4 rk4;

        auto f0 = F.eval(y, t);

        for (int init = 0; init < 3 && t < tEnd - 1e-12; ++init)
        {
            double h = std::min(dt, tEnd - t);
            y = rk4.step(F, y, t, h);
            t += h;
            res.t.push_back(t);
            res.y.push_back(y);
        }

        auto f1 = F.eval(res.y[1], res.t[1]);
        auto f2 = F.eval(res.y[2], res.t[2]);
        auto f3 = F.eval(res.y[3], res.t[3]);
        y = res.y[3];
        t = res.t[3];

        int s = 3;

        while (t < tEnd - 1e-12 * std::abs(dt) && s < maxSteps)
        {
            double h = std::min(dt, tEnd - t);

            auto yPred = linComb(y, h / 24.0,
                                 {{-9.0, &f0}, {37.0, &f1}, {-59.0, &f2}, {55.0, &f3}});
            auto fPred = F.eval(yPred, t + h);

            auto yNew = linComb(y, h / 24.0,
                                {{1.0, &f1}, {-5.0, &f2}, {19.0, &f3}, {9.0, &fPred}});

            f0 = f1;
            f1 = f2;
            f2 = f3;
            f3 = F.eval(yNew, t + h);
            y = yNew;
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