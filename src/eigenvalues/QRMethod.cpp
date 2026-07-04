// =============================================================================
// QRMethod.cpp
// Decomposição QR via Gram-Schmidt modificado + iteração A_{k+1} = R_k Q_k.
// =============================================================================
#include "eigenvalues/QRMethod.hpp"
#include "eigenvalues/HouseholderMethod.hpp"
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace nm
{

    // ---------------------------------------------------------------------------
    // Gram-Schmidt modificado: decompõe A = Q R
    // Q: ortogonal (colunas ortonormais), R: triangular superior
    // ---------------------------------------------------------------------------
    void QRMethod::qrDecompose(const Matrix &A, Matrix &Q, Matrix &R) const
    {
        int n = A.rows();
        Q = A; // Começa com as colunas de A; serão ortonormalizadas
        R = Matrix(n, n, 0.0);

        for (int j = 0; j < n; ++j)
        {
            // Extrai coluna j de Q
            std::vector<double> qj(n);
            for (int i = 0; i < n; ++i)
                qj[i] = Q(i, j);

            // Ortogonaliza em relação às colunas anteriores (Gram-Schmidt modificado)
            for (int k = 0; k < j; ++k)
            {
                std::vector<double> qk(n);
                for (int i = 0; i < n; ++i)
                    qk[i] = Q(i, k);

                double rjk = Matrix::dot(qk, qj); // r_{kj} = qk · qj
                R(k, j) = rjk;
                for (int i = 0; i < n; ++i)
                    qj[i] -= rjk * qk[i]; // qj ← qj - r_{kj} qk
            }

            // Normaliza qj → coluna j de Q; r_{jj} = ||qj||
            double norm = Matrix::vecNorm(qj);
            if (norm < 1e-14)
            {
                // Coluna linearmente dependente: insere vetor canônico
                for (int i = 0; i < n; ++i)
                    qj[i] = (i == j) ? 1.0 : 0.0;
                norm = 1.0;
            }
            R(j, j) = norm;
            for (int i = 0; i < n; ++i)
                Q(i, j) = qj[i] / norm;
        }
    }

    // ---------------------------------------------------------------------------
    // Convergência: máximo dos elementos subdiagonais < tol
    // ---------------------------------------------------------------------------
    bool QRMethod::isConverged(const Matrix &A, double tol) const
    {
        int n = A.rows();
        for (int i = 1; i < n; ++i)
            if (std::abs(A(i, i - 1)) > tol)
                return false;
        return true;
    }

    // ---------------------------------------------------------------------------
    // Iteração QR principal
    // ---------------------------------------------------------------------------
    QRResult QRMethod::solve(const Matrix &A, double tol, int maxIter) const
    {
        if (!A.isSquare())
            throw std::invalid_argument("QRMethod: matriz deve ser quadrada.");
        int n = A.rows();

        Matrix Ak = A;
        Matrix V = Matrix::identity(n); // Acumula autovetores: V = Q1 Q2 ...

        // Pré-tridiagonalização opcional (acelera convergência)
        if (useHouseholder_ && n > 2)
        {
            HouseholderMethod hh;
            auto res = hh.tridiagonalize(A);
            Ak = res.T;
            V = res.Q;
        }

        int iter = 0;
        bool converged = false;

        for (iter = 1; iter <= maxIter; ++iter)
        {
            Matrix Q, R;
            qrDecompose(Ak, Q, R);

            // ==============================
            // PRIMEIRA ITERAÇÃO
            // ==============================
            if (iter == 1)
            {
                std::cout << "\n=== Primeira iteração QR ===\n";

                std::cout << "\nMatriz Q1:\n";
                Q.print(); // ou seu método de impressão

                std::cout << "\nMatriz R1:\n";
                R.print();

                Matrix T1 = R * Q;

                std::cout << "\nT(1) = R1 * Q1:\n";
                T1.print();
            }

            // Atualização padrão
            Ak = R * Q;
            V = V * Q;

            if (isConverged(Ak, tol))
            {
                converged = true;
                break;
            }
        }

        // Diagonal de Ak = autovalores
        std::vector<double> eigenvalues(n);
        for (int i = 0; i < n; ++i)
            eigenvalues[i] = Ak(i, i);

        return {eigenvalues, V, iter, converged};
    }

} // namespace nm