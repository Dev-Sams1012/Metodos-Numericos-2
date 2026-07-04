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

            return ((-15.0 / 2.0) * f(x, 0, h) + 24.0 * f(x, 1, h) - 30.0 * f(x, 2, h) + 16.0 * f(x, 3, h) - (3.0 / 2.0) * f(x, 4, h)) / h3;

        case ThirdDerivMethod::CENTRAL_O4:

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

}