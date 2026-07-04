#include "integration/SingularityStrategies.hpp"
#include <cmath>
#include <stdexcept>

namespace nm
{

    static constexpr double PI = 3.14159265358979323846;

    SingularityIntegrator::SingularityIntegrator(std::shared_ptr<FunctionParser> f,
                                                 SingularityMethod method,
                                                 double T, double hStep)
        : Integrator(std::move(f)), method_(method), T_(T), hStep_(hStep) {}

    double SingularityIntegrator::compute(double a, double b) const
    {
        switch (method_)
        {
        case SingularityMethod::SINGLE_EXP:
            return singleExp(a, b);
        case SingularityMethod::DOUBLE_EXP:
            return doubleExp(a, b);
        }
        throw std::runtime_error("SingularityIntegrator: método desconhecido.");
    }

    double SingularityIntegrator::singleExp(double a, double b) const
    {
        const double mid = (a + b) / 2.0;
        const double half = (b - a) / 2.0;
        double result = 0.0;

        int N = static_cast<int>(2.0 * T_ / hStep_) + 1;
        for (int k = 0; k <= N; ++k)
        {
            double s = -T_ + k * hStep_;
            double ts = std::tanh(s);
            double sec2 = 1.0 - ts * ts;
            double xs = mid + half * ts;
            double dxds = half * sec2;

            double weight = (k == 0 || k == N) ? 0.5 : 1.0;
            result += weight * eval(xs) * dxds;
        }
        return result * hStep_;
    }

    double SingularityIntegrator::doubleExp(double a, double b) const
    {
        const double mid = (a + b) / 2.0;
        const double half = (b - a) / 2.0;
        double result = 0.0;

        int N = static_cast<int>(2.0 * T_ / hStep_) + 1;
        for (int k = 0; k <= N; ++k)
        {
            double s = -T_ + k * hStep_;
            double sinhs = std::sinh(s);
            double coshs = std::cosh(s);
            double arg = (PI / 2.0) * sinhs;
            double targ = std::tanh(arg);
            double carg = std::cosh(arg);
            double xs = mid + half * targ;

            double dxds = half * (PI / 2.0) * coshs / (carg * carg);

            double weight = (k == 0 || k == N) ? 0.5 : 1.0;
            result += weight * eval(xs) * dxds;
        }
        return result * hStep_;
    }

    std::string SingularityIntegrator::name() const
    {
        switch (method_)
        {
        case SingularityMethod::SINGLE_EXP:
            return "Singularidade — Exponencial Simples (tanh)";
        case SingularityMethod::DOUBLE_EXP:
            return "Singularidade — Exponencial Dupla (tanh-sinh)";
        }
        return "desconhecido";
    }

}