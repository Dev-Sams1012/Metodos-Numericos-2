#pragma once
// =============================================================================
// SVDSolver.hpp
// Decomposição em Valores Singulares: A = U Σ Vᵀ
//
// Estratégia via autovalores de AᵀA:
//   1. Calcula B = Aᵀ A  (simétrica e semi-definida positiva)
//   2. Usa Jacobi (ou QR) para obter B = V D Vᵀ
//      → valores singulares σᵢ = √(λᵢ),  λᵢ = diagonal de D
//   3. Calcula U: uᵢ = A vᵢ / σᵢ  (para σᵢ > ε_tol)
//
// Resultado: A = U Σ Vᵀ, com
//   U  ∈ ℝ^{m×r}, colunas ortonormais (vetores singulares esquerdos)
//   Σ  ∈ ℝ^{r×r}, diagonal com σ₁ ≥ σ₂ ≥ … ≥ σᵣ > 0
//   V  ∈ ℝ^{n×r}, colunas ortonormais (vetores singulares direitos)
//   r  = posto numérico (número de σᵢ > tol)
// =============================================================================
#include "Matrix.hpp"
#include <vector>
#include <string>

namespace nm {

struct SVDResult {
    Matrix              U;               ///< Vetores singulares esquerdos (m×r)
    std::vector<double> singularValues;  ///< σ₁ ≥ … ≥ σᵣ (ordenados)
    Matrix              V;               ///< Vetores singulares direitos  (n×r)
    int                 rank;            ///< Posto numérico
};

class SVDSolver {
public:
    /// @param tol Limiar para considerar valor singular como zero.
    explicit SVDSolver(double tol = 1e-10) : tol_(tol) {}

    /// @brief Calcula A = U Σ Vᵀ.
    /// @param A Matriz m×n qualquer (não precisa ser quadrada nem simétrica).
    SVDResult solve(const Matrix& A) const;

    std::string name() const {
        return "SVD — Decomposição em Valores Singulares (via autovalores de AᵀA)";
    }

private:
    double tol_;
};

} // namespace nm