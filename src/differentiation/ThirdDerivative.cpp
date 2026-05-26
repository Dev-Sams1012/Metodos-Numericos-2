// =============================================================================
// ThirdDerivative.cpp
// =============================================================================
#include "differentiation/ThirdDerivative.hpp"

namespace nm
{

    ThirdDerivative::ThirdDerivative(std::shared_ptr<FunctionParser> fp, ThirdDerivMethod method)
        : Differentiator(std::move(fp)), method_(method) {}

    double ThirdDerivative::compute(double x, double h) const
    {
        const double h3 = h * h * h;
        switch (method_)
        {
        case ThirdDerivMethod::FORWARD_O2:
            // f'''(xi) ≈ (1/h³) * [-15/2 fi + 24 f_{i+1} - 30 f_{i+2} + 16 f_{i+3} - 3/2 f_{i+4}]
            return ((-15.0 / 2.0) * f(x, 0, h) + 24.0 * f(x, 1, h) - 30.0 * f(x, 2, h) + 16.0 * f(x, 3, h) - (3.0 / 2.0) * f(x, 4, h)) / h3;

        case ThirdDerivMethod::CENTRAL_O4:
            // f'''(xi) ≈ (1 / 8h³) * [-f_{i+3} + 8f_{i+2} - 13f_{i+1} + 13f_{i-1} - 8f_{i-2} + f_{i-3}]
            return (-f(x, 3, h) + 8.0 * f(x, 2, h) - 13.0 * f(x, 1, h) + 13.0 * f(x, -1, h) - 8.0 * f(x, -2, h) + f(x, -3, h)) / (8.0 * h3);
        }
        throw std::runtime_error("ThirdDerivative: método desconhecido.");
    }

    std::string ThirdDerivative::name() const
    {
        switch (method_)
        {
        case ThirdDerivMethod::FORWARD_O2:
            return "f''' Forward  O(h²)";
        case ThirdDerivMethod::CENTRAL_O4:
            return "f''' Central  O(h⁴)";
        }
        return "desconhecido";
    }

} // namespace nm