// =============================================================================
// SecondDerivative.cpp
// =============================================================================
#include "differentiation/SecondDerivative.hpp"

namespace nm
{

    SecondDerivative::SecondDerivative(std::shared_ptr<FunctionParser> fp, SecondDerivMethod method)
        : Differentiator(std::move(fp)), method_(method) {}

    double SecondDerivative::compute(double x, double h) const
    {
        const double h2 = h * h;
        switch (method_)
        {
        case SecondDerivMethod::CENTRAL_3PT:
            // f''(xi) ≈ [f(xi+1) - 2f(xi) + f(xi-1)] / h²
            return (f(x, 1, h) - 2.0 * f(x, 0, h) + f(x, -1, h)) / h2;

        case SecondDerivMethod::CENTRAL_5PT:
            // f''(xi) ≈ [-f(xi+2) + 16f(xi+1) - 30f(xi) + 16f(xi-1) - f(xi-2)] / (12h²)
            return (-f(x, 2, h) + 16.0 * f(x, 1, h) - 30.0 * f(x, 0, h) + 16.0 * f(x, -1, h) - f(x, -2, h)) / (12.0 * h2);
        }
        throw std::runtime_error("SecondDerivative: método desconhecido.");
    }

    std::string SecondDerivative::name() const
    {
        switch (method_)
        {
        case SecondDerivMethod::CENTRAL_3PT:
            return "f'' Central 3pt O(h²)";
        case SecondDerivMethod::CENTRAL_5PT:
            return "f'' Central     O(h⁴)";
        }
        return "desconhecido";
    }

} // namespace nm