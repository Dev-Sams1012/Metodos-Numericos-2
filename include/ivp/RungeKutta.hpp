#pragma once
// =============================================================================
// RungeKutta.hpp
// Métodos de Runge-Kutta de passo simples explícitos.
//
//   RK2  (Heun / Trapézio):   O(Δt²)
//     k1 = F(Si, ti)
//     k2 = F(Si + Δt·k1, ti+Δt)
//     S_{i+1} = Si + (Δt/2)·(k1 + k2)
//
//   RK3  (Simpson 1/3, Kutta):  O(Δt³)
//     k1 = F(Si, ti)
//     k2 = F(Si + (Δt/2)·k1,        ti + Δt/2)
//     k3 = F(Si + Δt·(-k1 + 2·k2),  ti + Δt)
//     S_{i+1} = Si + (Δt/6)·(k1 + 4·k2 + k3)
//
//   RK4  (Clássico):  O(Δt⁴)
//     k1 = F(Si, ti)
//     k2 = F(Si + (Δt/2)·k1, ti + Δt/2)
//     k3 = F(Si + (Δt/2)·k2, ti + Δt/2)
//     k4 = F(Si + Δt·k3,     ti + Δt)
//     S_{i+1} = Si + (Δt/6)·(k1 + 2·k2 + 2·k3 + k4)
// =============================================================================
#include "ivp/IVPSolver.hpp"

namespace nm
{

// ---------------------------------------------------------------------------
class RungeKutta2 : public IVPSolver
{
public:
    RungeKutta2() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override { return "Runge-Kutta 2ª ordem  O(Δt²)  (Heun)"; }

    /// @brief Realiza um único passo RK2 — útil para inicialização de AB.
    std::vector<double> step(const SystemParser& F,
                              const std::vector<double>& y,
                              double t, double dt) const;
};

// ---------------------------------------------------------------------------
class RungeKutta3 : public IVPSolver
{
public:
    RungeKutta3() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override { return "Runge-Kutta 3ª ordem  O(Δt³)  (Kutta/Simpson)"; }

    std::vector<double> step(const SystemParser& F,
                              const std::vector<double>& y,
                              double t, double dt) const;
};

// ---------------------------------------------------------------------------
class RungeKutta4 : public IVPSolver
{
public:
    RungeKutta4() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override { return "Runge-Kutta 4ª ordem  O(Δt⁴)  (Clássico)"; }

    std::vector<double> step(const SystemParser& F,
                              const std::vector<double>& y,
                              double t, double dt) const;
};

} // namespace nm