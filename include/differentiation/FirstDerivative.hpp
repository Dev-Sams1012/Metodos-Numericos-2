#pragma once

#include "differentiation/Differentiator.hpp"

namespace nm
{

    enum class FirstDerivMethod
    {
        FORWARD_2PT,
        FORWARD_3PT,
        FORWARD_4PT,
        BACKWARD_2PT,
        BACKWARD_3PT,
        BACKWARD_4PT,
        CENTRAL_2PT,
        CENTRAL_4PT
    };

    class FirstDerivative : public Differentiator
    {
    public:
        FirstDerivative(std::shared_ptr<FunctionParser> f, FirstDerivMethod method);

        double compute(double x, double h) const override;
        std::string name() const override;

        void setMethod(FirstDerivMethod m) { method_ = m; }

    private:
        FirstDerivMethod method_;

        double forward2(double x, double h) const;
        double forward3(double x, double h) const;
        double forward4(double x, double h) const;
        double backward2(double x, double h) const;
        double backward3(double x, double h) const;
        double backward4(double x, double h) const;
        double central2(double x, double h) const;
        double central4(double x, double h) const;
    };

}