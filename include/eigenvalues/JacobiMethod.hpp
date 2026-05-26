#pragma once
// =============================================================================
// JacobiMethod.hpp
// Diagonalização iterativa de Jacobi para matrizes simétricas reais.
//
// A cada iteração:
//   1. Encontra o maior elemento fora da diagonal: a_{pq}
//   2. Calcula θ: tan(2θ) = -2 a_{pq} / (a_{pp} - a_{qq})
//      Se a_{pp} == a_{qq}: θ = π/4 (ou -π/4 se a_{pq} < 0)
//   3. Monta J_{pq} (rotação de Givens):
//      J_{pp} = J_{qq} = c = cos(θ)
//      J_{pq} = -s,  J_{qp} = s = sin(θ)
//   4. A' = Jᵀ A J
//   5. Acumula V ← V J  (colunas de V → autovetores)
//
// Convergência: ||fora diagonal||_F < ε
// =============================================================================
#include "eigenvalues/EigenSolver.hpp"
#include <vector>

namespace nm {

struct JacobiResult {
    std::vector<double> eigenvalues;   ///< Diagonal de A diagonalizada
    Matrix              eigenvectors;  ///< Colunas = autovetores
    int                 iterations;
    bool                converged;
};

class JacobiMethod {
public:
    JacobiMethod() = default;

    /// @brief Diagonaliza A e retorna todos os autovalores/autovetores.
    JacobiResult solve(const Matrix& A, double tol = 1e-10,
                       int maxIter = 10000) const;

    std::string name() const { return "Método de Jacobi"; }
};

} // namespace nm