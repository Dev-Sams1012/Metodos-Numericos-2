#pragma once

#include "integration/Integrator.hpp"
#include <cmath>

namespace nm
{

    enum class SingularityMethod
    {
        SINGLE_EXP,
        DOUBLE_EXP
    };

    class SingularityIntegrator : public Integrator
    {
    public:
        SingularityIntegrator(std::shared_ptr<FunctionParser> f,
                              SingularityMethod method,
                              double T = 3.5,
                              double hStep = 0.05);

        double compute(double a, double b) const override;
        std::string name() const override;

        void setT(double T) { T_ = T; }
        void setStep(double hStep) { hStep_ = hStep; }

    private:
        SingularityMethod method_;
        double T_;
        double hStep_;

        double singleExp(double a, double b) const;
        double doubleExp(double a, double b) const;
    };

}