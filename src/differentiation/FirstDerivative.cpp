// =============================================================================
// FirstDerivative.cpp
// =============================================================================
#include "differentiation/FirstDerivative.hpp"

namespace nm
{

    FirstDerivative::FirstDerivative(std::shared_ptr<FunctionParser> fp, FirstDerivMethod method)
        : Differentiator(std::move(fp)), method_(method) {}

    double FirstDerivative::compute(double x, double h) const
    {
        switch (method_)
        {
        case FirstDerivMethod::FORWARD_2PT:
            return forward2(x, h);
        case FirstDerivMethod::FORWARD_3PT:
            return forward3(x, h);
        case FirstDerivMethod::FORWARD_4PT:
            return forward4(x, h);
        case FirstDerivMethod::BACKWARD_2PT:
            return backward2(x, h);
        case FirstDerivMethod::BACKWARD_3PT:
            return backward3(x, h);
        case FirstDerivMethod::BACKWARD_4PT:
            return backward4(x, h);
        case FirstDerivMethod::CENTRAL_2PT:
            return central2(x, h);
        case FirstDerivMethod::CENTRAL_4PT:
            return central4(x, h);
        }
        throw std::runtime_error("FirstDerivative: método desconhecido.");
    }

    std::string FirstDerivative::name() const
    {
        switch (method_)
        {
        case FirstDerivMethod::FORWARD_2PT:
            return "f'  Forward  2pt  O(h)";
        case FirstDerivMethod::FORWARD_3PT:
            return "f'  Forward  3pt  O(h²)";
        case FirstDerivMethod::FORWARD_4PT:
            return "f'  Forward  4pt  O(h³)";
        case FirstDerivMethod::BACKWARD_2PT:
            return "f'  Backward 2pt  O(h)";
        case FirstDerivMethod::BACKWARD_3PT:
            return "f'  Backward 3pt  O(h²)";
        case FirstDerivMethod::BACKWARD_4PT:
            return "f'  Backward 4pt  O(h³)";
        case FirstDerivMethod::CENTRAL_2PT:
            return "f'  Central  2pt  O(h²)";
        case FirstDerivMethod::CENTRAL_4PT:
            return "f'  Central  4pt  O(h⁴)";
        }
        return "desconhecido";
    }

    // ---- Implementações das fórmulas -------------------------------------------

    // f'(xi) ≈ [f(xi+1) - f(xi)] / h
    double FirstDerivative::forward2(double x, double h) const
    {
        return (f(x, 1, h) - f(x, 0, h)) / h;
    }

    // f'(xi) ≈ [-f(xi+2) + 4f(xi+1) - 3f(xi)] / (2h)
    double FirstDerivative::forward3(double x, double h) const
    {
        return (-f(x, 2, h) + 4.0 * f(x, 1, h) - 3.0 * f(x, 0, h)) / (2.0 * h);
    }

    // f'(xi) ≈ [2f(xi+3) - 9f(xi+2) + 18f(xi+1) - 11f(xi)] / (6h)
    double FirstDerivative::forward4(double x, double h) const
    {
        return (2.0 * f(x, 3, h) - 9.0 * f(x, 2, h) + 18.0 * f(x, 1, h) - 11.0 * f(x, 0, h)) / (6.0 * h);
    }

    // f'(xi) ≈ [f(xi) - f(xi-1)] / h
    double FirstDerivative::backward2(double x, double h) const
    {
        return (f(x, 0, h) - f(x, -1, h)) / h;
    }

    // f'(xi) ≈ [3f(xi) - 4f(xi-1) + f(xi-2)] / (2h)
    double FirstDerivative::backward3(double x, double h) const
    {
        return (3.0 * f(x, 0, h) - 4.0 * f(x, -1, h) + f(x, -2, h)) / (2.0 * h);
    }

    // f'(xi) ≈ [11f(xi) - 18f(xi-1) + 9f(xi-2) - 2f(xi-3)] / (6h)
    double FirstDerivative::backward4(double x, double h) const
    {
        return (11.0 * f(x, 0, h) - 18.0 * f(x, -1, h) + 9.0 * f(x, -2, h) - 2.0 * f(x, -3, h)) / (6.0 * h);
    }

    // f'(xi) ≈ [f(xi+1) - f(xi-1)] / (2h)
    double FirstDerivative::central2(double x, double h) const
    {
        return (f(x, 1, h) - f(x, -1, h)) / (2.0 * h);
    }

    // f'(xi) ≈ [-f(xi+2) + 8f(xi+1) - 8f(xi-1) + f(xi-2)] / (12h)
    double FirstDerivative::central4(double x, double h) const
    {
        return (-f(x, 2, h) + 8.0 * f(x, 1, h) - 8.0 * f(x, -1, h) + f(x, -2, h)) / (12.0 * h);
    }

} // namespace nm