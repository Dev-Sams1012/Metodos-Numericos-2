#pragma once

#include "integration/Integrator.hpp"
#include <vector>

namespace nm
{

    class GaussLegendre : public Integrator
    {
    public:
        GaussLegendre(std::shared_ptr<FunctionParser> f, int n = 5);

        double compute(double a, double b) const override;
        std::string name() const override;

        void setPoints(int n);
        int points() const { return n_; }

    private:
        int n_;
        std::vector<double> nodes_;
        std::vector<double> weights_;

        void loadTable(int n);
    };

}