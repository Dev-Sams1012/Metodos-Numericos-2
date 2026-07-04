#pragma once

#include "integration/Integrator.hpp"

namespace nm
{

    enum class NewtonCotesMethod
    {

        TRAPEZOID,
        SIMPSON_1_3,
        SIMPSON_3_8,

        OPEN_TRAPEZOID,
        MILNE
    };

    class NewtonCotes : public Integrator
    {
    public:
        NewtonCotes(std::shared_ptr<FunctionParser> f, NewtonCotesMethod method);

        double compute(double a, double b) const override;
        std::string name() const override;

        void setMethod(NewtonCotesMethod m) { method_ = m; }

    private:
        NewtonCotesMethod method_;
    };

}