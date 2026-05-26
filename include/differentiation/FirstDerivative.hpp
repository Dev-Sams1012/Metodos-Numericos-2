#pragma once
// =============================================================================
// FirstDerivative.hpp
// Aproximações de diferenças finitas para a primeira derivada f'(xi).
//
// Métodos implementados:
//   Forward  2 pt: [f(xi+1) - f(xi)] / h
//   Forward  3 pt: [-f(xi+2) + 4f(xi+1) - 3f(xi)] / (2h)
//   Forward  4 pt: [2f(xi+3) - 9f(xi+2) + 18f(xi+1) - 11f(xi)] / (6h)
//   Backward 2 pt: [f(xi) - f(xi-1)] / h
//   Backward 3 pt: [3f(xi) - 4f(xi-1) + f(xi-2)] / (2h)
//   Backward 4 pt: [11f(xi) - 18f(xi-1) + 9f(xi-2) - 2f(xi-3)] / (6h)
//   Central  2 pt: [f(xi+1) - f(xi-1)] / (2h)
//   Central  4 pt: [-f(xi+2) + 8f(xi+1) - 8f(xi-1) + f(xi-2)] / (12h)
// =============================================================================
#include "differentiation/Differentiator.hpp"

namespace nm
{

    /// @brief Enumeração de todos os esquemas para f'(x).
    enum class FirstDerivMethod
    {
        FORWARD_2PT,  ///< Forward 2 pontos  O(h)
        FORWARD_3PT,  ///< Forward 3 pontos  O(h²)
        FORWARD_4PT,  ///< Forward 4 pontos  O(h³)
        BACKWARD_2PT, ///< Backward 2 pontos O(h)
        BACKWARD_3PT, ///< Backward 3 pontos O(h²)
        BACKWARD_4PT, ///< Backward 4 pontos O(h³)
        CENTRAL_2PT,  ///< Central  2 pontos O(h²)
        CENTRAL_4PT   ///< Central  4 pontos O(h⁴)
    };

    /// @brief Calcula f'(x) com o método escolhido via polimorfismo de seleção.
    class FirstDerivative : public Differentiator
    {
    public:
        FirstDerivative(std::shared_ptr<FunctionParser> f, FirstDerivMethod method);

        double compute(double x, double h) const override;
        std::string name() const override;

        void setMethod(FirstDerivMethod m) { method_ = m; }

    private:
        FirstDerivMethod method_;

        // Implementações internas de cada fórmula
        double forward2(double x, double h) const;
        double forward3(double x, double h) const;
        double forward4(double x, double h) const;
        double backward2(double x, double h) const;
        double backward3(double x, double h) const;
        double backward4(double x, double h) const;
        double central2(double x, double h) const;
        double central4(double x, double h) const;
    };

} // namespace nm