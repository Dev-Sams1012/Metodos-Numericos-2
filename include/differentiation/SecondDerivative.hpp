#pragma once
// =============================================================================
// SecondDerivative.hpp
// Aproximações de diferenças finitas para a segunda derivada f''(xi).
//
// Métodos implementados:
//   Central 3pt: [f(xi+1) - 2f(xi) + f(xi-1)] / h²
//   Central 5pt: [-f(xi+2) + 16f(xi+1) - 30f(xi) + 16f(xi-1) - f(xi-2)] / (12h²)
// =============================================================================
#include "differentiation/Differentiator.hpp"

namespace nm
{

    enum class SecondDerivMethod
    {
        CENTRAL_3PT, ///< O(h²)
        CENTRAL_5PT  ///< O(h⁴)
    };

    class SecondDerivative : public Differentiator
    {
    public:
        SecondDerivative(std::shared_ptr<FunctionParser> f, SecondDerivMethod method);

        double compute(double x, double h) const override;
        std::string name() const override;

        void setMethod(SecondDerivMethod m) { method_ = m; }

    private:
        SecondDerivMethod method_;

        double central3(double x, double h) const;
        double central5(double x, double h) const;
    };

} // namespace nm