#pragma once
// =============================================================================
// AdamsBashforth.hpp
// Métodos de passos múltiplos do tipo preditor-corretor (Adams-Bashforth-Moulton).
//
//   AB2  (2ª ordem):  O(Δt²)
//     Inicialização: 1 passo de RK2
//     Preditor:   S* = Si + (Δt/2)·(-Fi-1 + 3·Fi)
//     Corretor:   S  = Si + (Δt/2)·(Fi + Fi+1*)
//
//   AB3  (3ª ordem):  O(Δt³)
//     Inicialização: 2 passos de RK3
//     Preditor:   S* = Si + (Δt/12)·(5·Fi-2 - 16·Fi-1 + 23·Fi)
//     Corretor:   S  = Si + (Δt/12)·(-Fi-1 + 8·Fi + 5·Fi+1*)
//
//   AB4  (4ª ordem):  O(Δt⁴)
//     Inicialização: 3 passos de RK4
//     Preditor:   S* = Si + (Δt/24)·(-9·Fi-3 + 37·Fi-2 - 59·Fi-1 + 55·Fi)
//     Corretor:   S  = Si + (Δt/24)·(Fi-2 - 5·Fi-1 + 19·Fi + 9·Fi+1*)
// =============================================================================
#include "ivp/IVPSolver.hpp"

namespace nm
{

// ---------------------------------------------------------------------------
class AdamsBashforth2 : public IVPSolver
{
public:
    AdamsBashforth2() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override
    { return "Adams-Bashforth-Moulton 2ª ordem  O(Δt²)  (PC, init: RK2)"; }
};

// ---------------------------------------------------------------------------
class AdamsBashforth3 : public IVPSolver
{
public:
    AdamsBashforth3() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override
    { return "Adams-Bashforth-Moulton 3ª ordem  O(Δt³)  (PC, init: RK3)"; }
};

// ---------------------------------------------------------------------------
class AdamsBashforth4 : public IVPSolver
{
public:
    AdamsBashforth4() = default;

    IVPResult solve(const SystemParser& F,
                    const std::vector<double>& y0,
                    double t0, double tEnd, double dt,
                    int maxSteps = 200000) const override;

    std::string name() const override
    { return "Adams-Bashforth-Moulton 4ª ordem  O(Δt⁴)  (PC, init: RK4)"; }
};

} // namespace nm