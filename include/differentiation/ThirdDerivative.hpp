#pragma once

#include "differentiation/Differentiator.hpp"

namespace nm
{

    enum class ThirdDerivMethod
    {
        FORWARD_O2,
        CENTRAL_O4
    };

    class ThirdDerivative : public Differentiator
    {
    public:
        ThirdDerivative(std::shared_ptr<FunctionParser> f, ThirdDerivMethod method);

        double compute(double x, double h) const override;
        std::string name() const override;

        void setMethod(ThirdDerivMethod m) { method_ = m; }

    private:
        ThirdDerivMethod method_;

        double forwardO2(double x, double h) const;
        double centralO4(double x, double h) const;
    };

} 