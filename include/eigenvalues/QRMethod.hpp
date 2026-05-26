#pragma once
// =============================================================================
// QRMethod.hpp
// Método QR iterativo para todos os autovalores de uma matriz simétrica.
//
// Algoritmo:
//   1. (Opcional, recomendado) Tridiagonaliza A via Householder → T
//   2. Loop: A_k = Q_k R_k  (decomposição QR via Gram-Schmidt modificado)
//            A_{k+1} = R_k Q_k
//   3. Convergência quando os elementos subdiagonais de A ficam < ε
//
// Quando A já é tridiagonal, cada iteração QR tem custo O(n) em vez de O(n²).
// A diagonal de A convergida contém os autovalores.
// =============================================================================
#include "eigenvalues/EigenSolver.hpp"
#include <vector>

namespace nm {

struct QRResult {
    std::vector<double> eigenvalues;   ///< Diagonal da matriz convergida
    Matrix              eigenvectors;  ///< Colunas = autovetores acumulados
    int                 iterations;
    bool                converged;
};

class QRMethod {
public:
    /// @param useHouseholder Pré-tridiagonaliza com Householder antes do QR
    explicit QRMethod(bool useHouseholder = true)
        : useHouseholder_(useHouseholder) {}

    /// @brief Calcula todos os autovalores (e autovetores) de A.
    QRResult solve(const Matrix& A, double tol = 1e-10,
                   int maxIter = 10000) const;

    std::string name() const { return "Método QR (decomposição iterativa)"; }

private:
    bool useHouseholder_;

    /// @brief Decomposição QR via Gram-Schmidt modificado.
    /// @param A  Matriz de entrada (n×n).
    /// @param Q  Saída: matriz ortogonal.
    /// @param R  Saída: matriz triangular superior.
    void qrDecompose(const Matrix& A, Matrix& Q, Matrix& R) const;

    /// @brief Verifica convergência: ||subdiagonal||_∞ < tol
    bool isConverged(const Matrix& A, double tol) const;
};

} // namespace nm