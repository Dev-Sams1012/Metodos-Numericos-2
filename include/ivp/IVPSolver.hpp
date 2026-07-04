#pragma once

#include "core/SystemParser.hpp"
#include <vector>
#include <string>

namespace nm
{

    struct IVPResult
    {
        std::vector<double> t;
        std::vector<std::vector<double>> y;
        int steps;
        bool completed;
    };

    class IVPSolver
    {
    public:
        virtual ~IVPSolver() = default;

        IVPSolver(const IVPSolver &) = delete;
        IVPSolver &operator=(const IVPSolver &) = delete;

        virtual IVPResult solve(const SystemParser &F,
                                const std::vector<double> &y0,
                                double t0, double tEnd, double dt,
                                int maxSteps = 200000) const = 0;

        virtual std::string name() const = 0;

    protected:
        IVPSolver() = default;

        static std::vector<double> axpy(const std::vector<double> &a,
                                        double alpha,
                                        const std::vector<double> &b)
        {
            int n = (int)a.size();
            std::vector<double> r(n);
            for (int i = 0; i < n; ++i)
                r[i] = a[i] + alpha * b[i];
            return r;
        }

        static double norm(const std::vector<double> &v)
        {
            double s = 0.0;
            for (double x : v)
                s += x * x;
            return std::sqrt(s);
        }

        static double relDiff(const std::vector<double> &a, const std::vector<double> &b)
        {
            int n = (int)a.size();
            double diff = 0.0;
            for (int i = 0; i < n; ++i)
                diff += (a[i] - b[i]) * (a[i] - b[i]);
            return std::sqrt(diff) / (1.0 + norm(b));
        }
    };

}