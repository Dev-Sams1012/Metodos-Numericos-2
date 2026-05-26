// =============================================================================
// SVDMethod.cpp
// =============================================================================
#include "eigenvalues/SVDMethod.hpp"
#include "eigenvalues/JacobiMethod.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace nm {

SVDResult SVDSolver::solve(const Matrix& A) const {
    int m = A.rows();
    int n = A.cols();

    // ------------------------------------------------------------------
    // Passo 1: B = Aᵀ A  (n×n, simétrica semidefinida positiva)
    // ------------------------------------------------------------------
    Matrix At = A.transpose();  // n×m
    Matrix B  = At * A;         // n×n

    // ------------------------------------------------------------------
    // Passo 2: Diagonaliza B via Jacobi → B = V D Vᵀ
    //   λᵢ = autovalores de B,  V = autovetores
    // ------------------------------------------------------------------
    JacobiMethod jacobi;
    JacobiResult jr = jacobi.solve(B, tol_ * tol_, 50000);

    int nEig = (int)jr.eigenvalues.size();

    // ------------------------------------------------------------------
    // Passo 3: Ordena σᵢ = √(λᵢ) de forma decrescente
    // ------------------------------------------------------------------
    // Índices de ordenação por autovalor descendente
    std::vector<int> idx(nEig);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return jr.eigenvalues[a] > jr.eigenvalues[b];
    });

    // Monta σ e V ordenados; conta posto
    std::vector<double> sigmas;
    sigmas.reserve(nEig);
    int rank = 0;
    for (int k : idx) {
        double lam = jr.eigenvalues[k];
        double sig = (lam > tol_ * tol_) ? std::sqrt(std::max(0.0, lam)) : 0.0;
        sigmas.push_back(sig);
        if (sig > tol_) ++rank;
    }

    // Monta V ordenada (n×r, só colunas com σ > tol)
    Matrix V(n, rank, 0.0);
    for (int col = 0; col < rank; ++col) {
        int srcCol = idx[col];
        for (int row = 0; row < n; ++row)
            V(row, col) = jr.eigenvectors(row, srcCol);
    }

    // ------------------------------------------------------------------
    // Passo 4: U = A V Σ⁻¹  (m×r)
    //   uᵢ = A vᵢ / σᵢ
    // ------------------------------------------------------------------
    Matrix U(m, rank, 0.0);
    for (int col = 0; col < rank; ++col) {
        double sig = sigmas[col];
        // Extrai vᵢ (coluna col de V)
        std::vector<double> vi(n);
        for (int row = 0; row < n; ++row) vi[row] = V(row, col);
        // uᵢ = A * vᵢ / σᵢ
        std::vector<double> ui = A.multiply(vi);
        for (int row = 0; row < m; ++row)
            U(row, col) = ui[row] / sig;
    }

    // Extrai apenas os r valores singulares não-nulos
    std::vector<double> svs(sigmas.begin(), sigmas.begin() + rank);

    return { U, svs, V, rank };
}

} // namespace nm