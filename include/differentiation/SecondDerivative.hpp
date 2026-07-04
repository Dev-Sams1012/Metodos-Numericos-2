#pragma once

#include "differentiation/Differentiator.hpp"

namespace nm
{

    enum class SecondDerivMethod
    {
        CENTRAL_3PT,
        CENTRAL_5PT
    };

    class SecondDerivative : public Differentiator
    {
    public:
        SecondDerivative(std::shared_ptr<FunctionParser> f, SecondDerivMethod method);

        double compute(double x, double h) const override;
        std::string name() const override;

        void setMethod(SecondDerivMethod m) { method_ = m; }

    private:
        SecondDerivMethod method_;

        double central3(double x, double h) const;
        double central5(double x, double h) const;
    };

} 