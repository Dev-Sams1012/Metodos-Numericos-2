#include "eigenvalues/QRMethod.hpp"
#include "eigenvalues/HouseholderMethod.hpp"
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace nm
{

    void QRMethod::qrDecompose(const Matrix &A, Matrix &Q, Matrix &R) const
    {
        int n = A.rows();
        Q = A;
        R = Matrix(n, n, 0.0);

        for (int j = 0; j < n; ++j)
        {

            std::vector<double> qj(n);
            for (int i = 0; i < n; ++i)
                qj[i] = Q(i, j);

            for (int k = 0; k < j; ++k)
            {
                std::vector<double> qk(n);
                for (int i = 0; i < n; ++i)
                    qk[i] = Q(i, k);

                double rjk = Matrix::dot(qk, qj);
                R(k, j) = rjk;
                for (int i = 0; i < n; ++i)
                    qj[i] -= rjk * qk[i];
            }

            double norm = Matrix::vecNorm(qj);
            if (norm < 1e-14)
            {

                for (int i = 0; i < n; ++i)
                    qj[i] = (i == j) ? 1.0 : 0.0;
                norm = 1.0;
            }
            R(j, j) = norm;
            for (int i = 0; i < n; ++i)
                Q(i, j) = qj[i] / norm;
        }
    }

    bool QRMethod::isConverged(const Matrix &A, double tol) const
    {
        int n = A.rows();
        for (int i = 1; i < n; ++i)
            if (std::abs(A(i, i - 1)) > tol)
                return false;
        return true;
    }

    QRResult QRMethod::solve(const Matrix &A, double tol, int maxIter) const
    {
        if (!A.isSquare())
            throw std::invalid_argument("QRMethod: matriz deve ser quadrada.");
        int n = A.rows();

        Matrix Ak = A;
        Matrix V = Matrix::identity(n);

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

            if (iter == 1)
            {
                std::cout << "\n=== Primeira iteração QR ===\n";

                std::cout << "\nMatriz Q1:\n";
                Q.print();

                std::cout << "\nMatriz R1:\n";
                R.print();

                Matrix T1 = R * Q;

                std::cout << "\nT(1) = R1 * Q1:\n";
                T1.print();
            }

            Ak = R * Q;
            V = V * Q;

            if (isConverged(Ak, tol))
            {
                converged = true;
                break;
            }
        }

        std::vector<double> eigenvalues(n);
        for (int i = 0; i < n; ++i)
            eigenvalues[i] = Ak(i, i);

        return {eigenvalues, V, iter, converged};
    }

}