#pragma once
// =============================================================================
// PowerMethod.hpp
// Três variantes do Método da Potência:
//
//   PowerRegular:   Autovalor DOMINANTE (|λ₁| máximo)
//     v_k = A * x_{k-1}
//     λ ≈ (x_{k-1}ᵀ v_k) / (x_{k-1}ᵀ x_{k-1})
//     x_k = v_k / ||v_k||
//
//   PowerInverse:   Autovalor de MENOR módulo
//     Aplica PowerRegular sobre A⁻¹ (resolve Ax = b via LU a cada passo)
//
//   PowerShifted:   Autovalor mais PRÓXIMO de μ
//     Aplica PowerInverse sobre (A - μI)
// =============================================================================
#include "eigenvalues/EigenSolver.hpp"

namespace nm {

// ---- Potência Regular ------------------------------------------------------
class PowerRegular : public EigenSolver {
public:
    PowerRegular() = default;
    EigenResult solve(const Matrix& A, double tol = 1e-10,
                      int maxIter = 1000) const override;
    std::string name() const override { return "Método da Potência Regular"; }
};

// ---- Potência Inverso -------------------------------------------------------
class PowerInverse : public EigenSolver {
public:
    PowerInverse() = default;
    EigenResult solve(const Matrix& A, double tol = 1e-10,
                      int maxIter = 1000) const override;
    std::string name() const override { return "Método da Potência Inverso"; }
};

// ---- Potência com Deslocamento ----------------------------------------------
class PowerShifted : public EigenSolver {
public:
    /// @param mu Estimativa do autovalor alvo (deslocamento).
    explicit PowerShifted(double mu = 0.0) : mu_(mu) {}

    void setShift(double mu) { mu_ = mu; }
    double shift() const     { return mu_; }

    EigenResult solve(const Matrix& A, double tol = 1e-10,
                      int maxIter = 1000) const override;
    std::string name() const override;

private:
    double mu_;
};

} // namespace nm