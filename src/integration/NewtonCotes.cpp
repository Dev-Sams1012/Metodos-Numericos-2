// =============================================================================
// NewtonCotes.cpp
// =============================================================================
#include "integration/NewtonCotes.hpp"

namespace nm
{

    NewtonCotes::NewtonCotes(std::shared_ptr<FunctionParser> f, NewtonCotesMethod method)
        : Integrator(std::move(f)), method_(method) {}

    double NewtonCotes::compute(double a, double b) const
    {
        switch (method_)
        {

            // ------ FECHADAS ---------------------------------------------------

        case NewtonCotesMethod::TRAPEZOID:
        {
            // (Δx/2)[f(x0) + f(x1)], Δx = b - a
            double dx = b - a;
            return (dx / 2.0) * (eval(a) + eval(b));
        }

        case NewtonCotesMethod::SIMPSON_1_3:
        {
            // h = (b-a)/2;  (h/3)[f(x0) + 4f(x1) + f(x2)]
            double h = (b - a) / 2.0;
            double x1 = a + h;
            return (h / 3.0) * (eval(a) + 4.0 * eval(x1) + eval(b));
        }

        case NewtonCotesMethod::SIMPSON_3_8:
        {
            // h = (b-a)/3;  (3h/8)[f(x0) + 3f(x1) + 3f(x2) + f(x3)]
            double h = (b - a) / 3.0;
            double x1 = a + h;
            double x2 = a + 2.0 * h;
            return (3.0 * h / 8.0) * (eval(a) + 3.0 * eval(x1) + 3.0 * eval(x2) + eval(b));
        }

            // ------ ABERTAS ----------------------------------------------------

        case NewtonCotesMethod::OPEN_TRAPEZOID:
        {
            // h = (b-a)/3;  (3h/2)[f(x1) + f(x2)]
            // x1 = a+h, x2 = a+2h  (pontos internos; extremos não usados)
            double h = (b - a) / 3.0;
            double x1 = a + h;
            double x2 = a + 2.0 * h;
            return (3.0 * h / 2.0) * (eval(x1) + eval(x2));
        }

        case NewtonCotesMethod::MILNE:
        {
            // h = (b-a)/4;  (4h/3)[2f(x1) - f(x2) + 2f(x3)]
            // x1=a+h, x2=a+2h, x3=a+3h
            double h = (b - a) / 4.0;
            double x1 = a + h;
            double x2 = a + 2.0 * h;
            double x3 = a + 3.0 * h;
            return (4.0 * h / 3.0) * (2.0 * eval(x1) - eval(x2) + 2.0 * eval(x3));
        }
        }
        throw std::runtime_error("NewtonCotes: método desconhecido.");
    }

    std::string NewtonCotes::name() const
    {
        switch (method_)
        {
        case NewtonCotesMethod::TRAPEZOID:
            return "Newton-Cotes — Trapézio (fechado)";
        case NewtonCotesMethod::SIMPSON_1_3:
            return "Newton-Cotes — Simpson 1/3";
        case NewtonCotesMethod::SIMPSON_3_8:
            return "Newton-Cotes — Simpson 3/8";
        case NewtonCotesMethod::OPEN_TRAPEZOID:
            return "Newton-Cotes — Trapézio Aberto";
        case NewtonCotesMethod::MILNE:
            return "Newton-Cotes — Regra de Milne";
        }
        return "desconhecido";
    }

} // namespace nm