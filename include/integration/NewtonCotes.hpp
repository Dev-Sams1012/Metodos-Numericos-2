#pragma once
// =============================================================================
// NewtonCotes.hpp
// Regras de Newton-Cotes (fechadas e abertas):
//
//   FECHADAS:
//     Trapézio:    (h/2) * [f(x0) + f(x1)],                     h = (b-a)
//     Simpson 1/3: (h/3) * [f(x0) + 4f(x1) + f(x2)],            h = (b-a)/2
//     Simpson 3/8: (3h/8) * [f(x0) + 3f(x1) + 3f(x2) + f(x3)],  h = (b-a)/3
//
//   ABERTAS:
//     Trapézio Aberto: (3h/2) * [f(x1) + f(x2)],                h = (b-a)/3
//     Milne (Aberta):  (4h/3) * [2f(x1) - f(x2) + 2f(x3)],      h = (b-a)/4
// =============================================================================
#include "integration/Integrator.hpp"

namespace nm
{

    enum class NewtonCotesMethod
    {
        // Fechadas
        TRAPEZOID,   ///< Regra do Trapézio
        SIMPSON_1_3, ///< Simpson 1/3
        SIMPSON_3_8, ///< Simpson 3/8

        // Abertas
        OPEN_TRAPEZOID, ///< Trapézio Aberto
        MILNE           ///< Regra de Milne (Aberta)
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

} // namespace nm