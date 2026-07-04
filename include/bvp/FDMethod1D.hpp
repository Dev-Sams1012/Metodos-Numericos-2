#pragma once

#include "bvp/BVPSolver.hpp"

namespace nm
{

    enum class BCType
    {
        DIRICHLET,
        NEUMANN
    };

    struct BoundaryCondition
    {
        BCType type = BCType::DIRICHLET;
        double value = 0.0;

        static BoundaryCondition dirichlet(double v)
        {
            return {BCType::DIRICHLET, v};
        }

        static BoundaryCondition neumann(double v)
        {
            return {BCType::NEUMANN, v};
        }
    };

    class FDMethod1D : public BVPSolver
    {
    public:
        explicit FDMethod1D(double alpha = 1.0,
                            double beta = 0.0,
                            double gamma = 0.0)
            : alpha_(alpha), beta_(beta), gamma_(gamma) {}

        BVPResult solve(const std::shared_ptr<FunctionParser> &f,
                        double a, double b,
                        const BoundaryCondition &bcA,
                        const BoundaryCondition &bcB,
                        int n) const;

        BVPResult solve(const std::shared_ptr<FunctionParser> &f,
                        double a, double b,
                        double ua, double ub,
                        int n) const
        {
            return solve(f, a, b,
                         BoundaryCondition::dirichlet(ua),
                         BoundaryCondition::dirichlet(ub),
                         n);
        }

        std::string name() const override
        {
            return "Diferenças Finitas 1D  (diferenças centrais O(h²), CC mistas D/N)";
        }

        void setAlpha(double v) { alpha_ = v; }
        void setBeta(double v) { beta_ = v; }
        void setGamma(double v) { gamma_ = v; }

    private:
        double alpha_;
        double beta_;
        double gamma_;
    };

}