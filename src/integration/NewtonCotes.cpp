#include "integration/NewtonCotes.hpp"

namespace nm
{

    NewtonCotes::NewtonCotes(std::shared_ptr<FunctionParser> f, NewtonCotesMethod method)
        : Integrator(std::move(f)), method_(method) {}

    double NewtonCotes::compute(double a, double b) const
    {
        switch (method_)
        {

        case NewtonCotesMethod::TRAPEZOID:
        {

            double dx = b - a;
            return (dx / 2.0) * (eval(a) + eval(b));
        }

        case NewtonCotesMethod::SIMPSON_1_3:
        {

            double h = (b - a) / 2.0;
            double x1 = a + h;
            return (h / 3.0) * (eval(a) + 4.0 * eval(x1) + eval(b));
        }

        case NewtonCotesMethod::SIMPSON_3_8:
        {

            double h = (b - a) / 3.0;
            double x1 = a + h;
            double x2 = a + 2.0 * h;
            return (3.0 * h / 8.0) * (eval(a) + 3.0 * eval(x1) + 3.0 * eval(x2) + eval(b));
        }

        case NewtonCotesMethod::OPEN_TRAPEZOID:
        {

            double h = (b - a) / 3.0;
            double x1 = a + h;
            double x2 = a + 2.0 * h;
            return (3.0 * h / 2.0) * (eval(x1) + eval(x2));
        }

        case NewtonCotesMethod::MILNE:
        {

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

}