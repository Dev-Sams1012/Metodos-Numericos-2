// =============================================================================
// PowerMethod.cpp
// =============================================================================
#include "eigenvalues/PowerMethod.hpp"
#include <cmath>
#include <numeric>
#include <sstream>
#include <random>

namespace nm {

// ---------- Helpers locais --------------------------------------------------
static std::vector<double> makeInitialVector(int n) {
    // Vetor inicial normalizado com componentes unitárias
    std::vector<double> v(n, 1.0);
    return Matrix::normalize(v);
}

// ===========================================================================
// PowerRegular — autovalor dominante
// ===========================================================================
EigenResult PowerRegular::solve(const Matrix& A, double tol, int maxIter) const {
    if (!A.isSquare())
        throw std::invalid_argument("PowerRegular: matriz deve ser quadrada.");
    int n = A.rows();

    std::vector<double> x = makeInitialVector(n);
    double lambda = 0.0;
    int iter = 0;
    bool converged = false;

    for (iter = 1; iter <= maxIter; ++iter) {
        // v_k = A * x_{k-1}
        std::vector<double> v = A.multiply(x);

        // λ ≈ (x_{k-1}ᵀ v_k) / (x_{k-1}ᵀ x_{k-1})
        double num = Matrix::dot(x, v);
        double den = Matrix::dot(x, x);
        double lambdaNew = num / den;

        // x_k = v_k / ||v_k||
        double normV = Matrix::vecNorm(v);
        if (normV < 1e-14) break; // Vetor nulo
        std::vector<double> xNew(n);
        for (int i = 0; i < n; ++i) xNew[i] = v[i] / normV;

        if (std::abs(lambdaNew - lambda) < tol) {
            lambda    = lambdaNew;
            x         = xNew;
            converged = true;
            break;
        }
        lambda = lambdaNew;
        x      = xNew;
    }

    return { lambda, x, iter, converged };
}

// ===========================================================================
// PowerInverse — menor autovalor (via A⁻¹)
// ===========================================================================
EigenResult PowerInverse::solve(const Matrix& A, double tol, int maxIter) const {
    if (!A.isSquare())
        throw std::invalid_argument("PowerInverse: matriz deve ser quadrada.");
    int n = A.rows();

    // Pré-fatoração LU de A (não muda a cada iteração → eficiente)
    Matrix LU = A; // cópia para decompor
    std::vector<int> piv;
    if (!LU.luDecompose(piv))
        throw std::runtime_error("PowerInverse: matriz singular, não há autovalor nulo.");

    std::vector<double> x = makeInitialVector(n);
    double lambda = 0.0;
    int iter = 0;
    bool converged = false;

    for (iter = 1; iter <= maxIter; ++iter) {
        // Resolve A * v = x  (equivale a multiplicar por A⁻¹)
        std::vector<double> v = LU.luSolve(x);

        // λ de A⁻¹ → autovalor de A = 1/λ_inv
        double num = Matrix::dot(x, v);
        double den = Matrix::dot(x, x);
        double lambdaInv = num / den;
        double lambdaNew = 1.0 / lambdaInv;

        double normV = Matrix::vecNorm(v);
        if (normV < 1e-14) break;
        std::vector<double> xNew(n);
        for (int i = 0; i < n; ++i) xNew[i] = v[i] / normV;

        if (std::abs(lambdaNew - lambda) < tol) {
            lambda    = lambdaNew;
            x         = xNew;
            converged = true;
            break;
        }
        lambda = lambdaNew;
        x      = xNew;
    }

    return { lambda, x, iter, converged };
}

// ===========================================================================
// PowerShifted — autovalor mais próximo de μ
// ===========================================================================
EigenResult PowerShifted::solve(const Matrix& A, double tol, int maxIter) const {
    if (!A.isSquare())
        throw std::invalid_argument("PowerShifted: matriz deve ser quadrada.");
    int n = A.rows();

    // B = A - μI
    Matrix B = A;
    for (int i = 0; i < n; ++i)
        B(i, i) -= mu_;

    // Pré-fatoração LU de B
    Matrix LU = B;
    std::vector<int> piv;
    if (!LU.luDecompose(piv))
        throw std::runtime_error("PowerShifted: (A - μI) singular; μ pode ser um autovalor exato.");

    std::vector<double> x = makeInitialVector(n);
    double lambda = 0.0;
    int iter = 0;
    bool converged = false;

    for (iter = 1; iter <= maxIter; ++iter) {
        // Resolve (A - μI) v = x
        std::vector<double> v = LU.luSolve(x);

        double num = Matrix::dot(x, v);
        double den = Matrix::dot(x, x);
        double lambdaShift = num / den;
        // Autovalor de A = 1/λ_shift + μ
        double lambdaNew = 1.0 / lambdaShift + mu_;

        double normV = Matrix::vecNorm(v);
        if (normV < 1e-14) break;
        std::vector<double> xNew(n);
        for (int i = 0; i < n; ++i) xNew[i] = v[i] / normV;

        if (std::abs(lambdaNew - lambda) < tol) {
            lambda    = lambdaNew;
            x         = xNew;
            converged = true;
            break;
        }
        lambda = lambdaNew;
        x      = xNew;
    }

    return { lambda, x, iter, converged };
}

std::string PowerShifted::name() const {
    std::ostringstream oss;
    oss << "Método da Potência com Deslocamento (μ=" << mu_ << ")";
    return oss.str();
}

} // namespace nm