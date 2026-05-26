#pragma once
// =============================================================================
// HouseholderMethod.hpp
// Transformação de Householder para reduzir uma matriz simétrica A a
// forma tridiagonal T = Hₙ₋₁ … H₂ H₁ A H₁ H₂ … Hₙ₋₁
//
// Reflexor de Householder: H = I - 2 w wᵀ / (wᵀ w)
// onde w = x ± ||x||e₁,  x = subcoluna abaixo da diagonal.
//
// A tridiagonalização é pré-condição para o método QR eficiente.
// =============================================================================
#include "Matrix.hpp"
#include <vector>
#include <string>

namespace nm {

struct HouseholderResult {
    Matrix              T;   ///< Matriz tridiagonal resultante
    Matrix              Q;   ///< Produto das transformações (Q tal que Q^T A Q = T)
    int                 steps;
};

class HouseholderMethod {
public:
    HouseholderMethod() = default;

    /// @brief Tridiagonaliza A via reflexores de Householder.
    HouseholderResult tridiagonalize(const Matrix& A) const;

    std::string name() const { return "Transformação de Householder (tridiagonalização)"; }
};

} // namespace nm