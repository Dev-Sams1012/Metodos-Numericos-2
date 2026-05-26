#pragma once
// =============================================================================
// SingularityStrategies.hpp
// Estratégias de mudança de variável para integrais com singularidades:
//
//   Exponencial Simples (tanh):
//     x(s) = (a+b)/2 + (b-a)/2 * tanh(s)
//     dx/ds = (b-a)/2 * sech²(s)
//     Domínio de s: (-∞, +∞) → truncado a [-T, T]
//
//   Exponencial Dupla (Double Exponential / tanh-sinh):
//     x(s) = (a+b)/2 + (b-a)/2 * tanh(π/2 * sinh(s))
//     dx/ds = (b-a)/2 * π/2 * cosh(s) / cosh²(π/2 * sinh(s))
//     Domínio de s: (-∞, +∞) → truncado a [-T, T]
//
// Ambas usam regra do trapézio uniforme em s com passo hStep.
// =============================================================================
#include "integration/Integrator.hpp"
#include <cmath>

namespace nm {

enum class SingularityMethod {
    SINGLE_EXP,  ///< Exponencial simples (tanh)
    DOUBLE_EXP   ///< Exponencial dupla (tanh-sinh)
};

class SingularityIntegrator : public Integrator {
public:
    /// @param f      Função a integrar.
    /// @param method Método de transformação.
    /// @param T      Truncamento do domínio de s (default 3.5).
    /// @param hStep  Passo em s (default 0.05).
    SingularityIntegrator(std::shared_ptr<FunctionParser> f,
                          SingularityMethod method,
                          double T     = 3.5,
                          double hStep = 0.05);

    double      compute(double a, double b) const override;
    std::string name()                      const override;

    void setT(double T)         { T_     = T;     }
    void setStep(double hStep)  { hStep_ = hStep; }

private:
    SingularityMethod method_;
    double            T_;      ///< Truncamento
    double            hStep_;  ///< Passo uniforme em s

    double singleExp(double a, double b) const;
    double doubleExp(double a, double b) const;
};

} // namespace nm