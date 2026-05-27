#pragma once
// =============================================================================
// EulerMethod.hpp
// Métodos de Euler para PVIs.
//
//   Euler Explícito:  S_{i+1} = S_i + Δt · F(S_i, t_i)           O(Δt)
//
//   Euler Implícito:  S_{i+1} = S_i + Δt · F(S_{i+1}, t_{i+1})   O(Δt)
//     Resolvido por iteração de ponto fixo:
//       S^(0)   = S_i + Δt·F(S_i, t_i)   (chute via Euler Explícito)
//       S^(k+1) = S_i + Δt·F(S^(k), t_{i+1})
//     até convergência.
// =============================================================================
#include "ivp/IVPSolver.hpp"

namespace nm
{

// ---------------------------------------------------------------------------
class EulerExplicit : public IVPSolver
{
public:
    EulerExplicit() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override { return "Euler Explícito  O(Δt)"; }
};

// ---------------------------------------------------------------------------
class EulerImplicit : public IVPSolver
{
public:
    /// @param tol     Tolerância de convergência do ponto fixo.
    /// @param maxIter Máximo de iterações do ponto fixo por passo.
    explicit EulerImplicit(double tol = 1e-10, int maxIter = 100)
        : tol_(tol), maxIter_(maxIter) {}

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override { return "Euler Implícito  O(Δt)  (ponto fixo)"; }

private:
    double tol_;
    int    maxIter_;
};

} // namespace nm