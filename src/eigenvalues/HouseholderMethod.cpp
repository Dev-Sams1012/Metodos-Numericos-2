// =============================================================================
// HouseholderMethod.cpp
// =============================================================================
#include "eigenvalues/HouseholderMethod.hpp"
#include <cmath>
#include <stdexcept>

namespace nm {

HouseholderResult HouseholderMethod::tridiagonalize(const Matrix& A) const {
    if (!A.isSquare())
        throw std::invalid_argument("Householder: matriz deve ser quadrada.");
    int n = A.rows();

    Matrix T = A;           // Cópia de trabalho — será tridiagonalizada
    Matrix Q = Matrix::identity(n); // Acumula transformações

    for (int k = 0; k < n - 2; ++k) {
        // Extrai subvetor x = T(k+1:n-1, k)  (coluna k abaixo da diagonal)
        int len = n - k - 1;
        std::vector<double> x(len);
        for (int i = 0; i < len; ++i)
            x[i] = T(k + 1 + i, k);

        // Calcula ||x||
        double normX = Matrix::vecNorm(x);
        if (normX < 1e-14) continue; // Subvetor já nulo

        // Escolhe sinal para minimizar cancelamento numérico
        double alpha = (x[0] >= 0.0) ? -normX : normX;

        // w = x - alpha * e₁
        std::vector<double> w(len);
        w[0] = x[0] - alpha;
        for (int i = 1; i < len; ++i) w[i] = x[i];

        // H = I - 2 w wᵀ / (wᵀ w)
        double wTw = Matrix::dot(w, w);
        if (wTw < 1e-14) continue;

        // Aplica H à submatriz: T' = H T H (eficientemente)
        // 1) T(k+1:, k:) = T - 2/wTw * w * (wᵀ T)
        // 2) T(:, k+1:)  = T - 2/wTw * (T w) * wᵀ
        // Fazemos de forma explícita (n pequeno para problemas típicos)

        // --- Aplica H pela esquerda: T ← H T
        for (int j = k; j < n; ++j) {
            // p_j = Σ_{i} w[i] * T(k+1+i, j)
            double p = 0.0;
            for (int i = 0; i < len; ++i)
                p += w[i] * T(k + 1 + i, j);
            p *= 2.0 / wTw;
            for (int i = 0; i < len; ++i)
                T(k + 1 + i, j) -= p * w[i];
        }

        // --- Aplica H pela direita: T ← T H
        for (int i = 0; i < n; ++i) {
            // q_i = Σ_{j} T(i, k+1+j) * w[j]
            double q = 0.0;
            for (int j = 0; j < len; ++j)
                q += T(i, k + 1 + j) * w[j];
            q *= 2.0 / wTw;
            for (int j = 0; j < len; ++j)
                T(i, k + 1 + j) -= q * w[j];
        }

        // --- Acumula em Q: Q ← Q H
        for (int i = 0; i < n; ++i) {
            double r = 0.0;
            for (int j = 0; j < len; ++j)
                r += Q(i, k + 1 + j) * w[j];
            r *= 2.0 / wTw;
            for (int j = 0; j < len; ++j)
                Q(i, k + 1 + j) -= r * w[j];
        }
    }

    return { T, Q, n - 2 };
}

} // namespace nm