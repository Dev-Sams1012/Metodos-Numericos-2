#pragma once
// =============================================================================
// FDMethod1D.hpp
// Método das Diferenças Finitas para PVCs unidimensionais.
//
// Problema resolvido (formulação geral):
//   α · d²u/dx² + β · du/dx + γ · u = f(x),   x ∈ [a, b]
//
// Condições de contorno suportadas (uma em cada extremidade):
//
//   Dirichlet : u(xCC)  = valor          (valor prescrito)
//   Neumann   : u'(xCC) = valor          (derivada prescrita)
//
// Todos os quatro combinações são válidas:
//   D–D   Dirichlet esq. + Dirichlet dir.
//   D–N   Dirichlet esq. + Neumann    dir.  ← exemplo da questão 5
//   N–D   Neumann    esq. + Dirichlet dir.
//   N–N   Neumann    esq. + Neumann    dir.
//
// Discretização no nó interior i (diferenças centrais, O(h²)):
//   u''(xi) ≈ [u_{i-1} - 2u_i + u_{i+1}] / h²
//   u' (xi) ≈ [u_{i+1} - u_{i-1}]        / (2h)
//
//   Equação discreta:
//     cL · u_{i-1}  +  cC · u_i  +  cR · u_{i+1} = f(xi)
//     cL = α/h² − β/(2h)
//     cC = −2α/h² + γ
//     cR = α/h² + β/(2h)
//
// Tratamento da condição de Neumann na borda direita (x = b = x_N):
//   Diferença retroativa O(h):   u'(b) ≈ (u_N − u_{N-1}) / h
//     → u_N = u_{N-1} + h · duB
//   Substituindo na equação do nó N-1 (último nó interior):
//     cL·u_{N-2} + cC·u_{N-1} + cR·(u_{N-1} + h·duB) = f(x_{N-1})
//     cL·u_{N-2} + (cC + cR)·u_{N-1} = f(x_{N-1}) − cR·h·duB
//
// Tratamento da condição de Neumann na borda esquerda (x = a = x_0):
//   Diferença progressiva O(h):  u'(a) ≈ (u_1 − u_0) / h
//     → u_0 = u_1 − h · duA
//   Substituindo na equação do nó 1 (primeiro nó interior):
//     cL·(u_1 − h·duA) + cC·u_1 + cR·u_2 = f(x_1)
//     (cL + cC)·u_1 + cR·u_2 = f(x_1) + cL·h·duA
//
// Sistema resultante: N_inc × N_inc, onde N_inc depende das CCs:
//   D–D → N-1 incógnitas (nós 1..N-1)
//   D–N → N   incógnitas (nós 1..N)    ← u_N passa a ser incógnita
//   N–D → N   incógnitas (nós 0..N-1)  ← u_0 passa a ser incógnita
//   N–N → N+1 incógnitas (nós 0..N)
//
// Passo 4: sistema resolvido via LU (Matrix::solve).
// =============================================================================
#include "bvp/BVPSolver.hpp"

namespace nm {

// ---------------------------------------------------------------------------
/// @brief Tipo de condição de contorno em uma extremidade.
enum class BCType {
    DIRICHLET,  ///< Valor prescrito:   u(x) = valor
    NEUMANN     ///< Derivada prescrita: u'(x) = valor
};

/// @brief Condição de contorno em uma extremidade.
struct BoundaryCondition {
    BCType type  = BCType::DIRICHLET;
    double value = 0.0;

    /// @brief Fábrica Dirichlet.
    static BoundaryCondition dirichlet(double v) {
        return { BCType::DIRICHLET, v };
    }
    /// @brief Fábrica Neumann.
    static BoundaryCondition neumann(double v) {
        return { BCType::NEUMANN, v };
    }
};
// ---------------------------------------------------------------------------

class FDMethod1D : public BVPSolver {
public:
    /// @brief Coeficientes do operador diferencial α·u'' + β·u' + γ·u = f(x).
    explicit FDMethod1D(double alpha = 1.0,
                        double beta  = 0.0,
                        double gamma = 0.0)
        : alpha_(alpha), beta_(beta), gamma_(gamma) {}

    // -----------------------------------------------------------------------
    /// @brief Resolve com condições de contorno genéricas (D ou N em cada lado).
    ///
    /// @param f    Função f(x). nullptr → f ≡ 0.
    /// @param a    Extremo esquerdo.
    /// @param b    Extremo direito.
    /// @param bcA  Condição em x = a  (Dirichlet ou Neumann).
    /// @param bcB  Condição em x = b  (Dirichlet ou Neumann).
    /// @param n    Número de divisões (h = (b-a)/n).
    BVPResult solve(const std::shared_ptr<FunctionParser>& f,
                    double a, double b,
                    const BoundaryCondition& bcA,
                    const BoundaryCondition& bcB,
                    int n) const;

    // -----------------------------------------------------------------------
    /// @brief Sobrecarga de conveniência: ambas as bordas Dirichlet.
    ///        Mantém compatibilidade com a interface anterior.
    BVPResult solve(const std::shared_ptr<FunctionParser>& f,
                    double a, double b,
                    double ua, double ub,
                    int n) const
    {
        return solve(f, a, b,
                     BoundaryCondition::dirichlet(ua),
                     BoundaryCondition::dirichlet(ub),
                     n);
    }

    std::string name() const override {
        return "Diferenças Finitas 1D  (diferenças centrais O(h²), CC mistas D/N)";
    }

    void setAlpha(double v) { alpha_ = v; }
    void setBeta (double v) { beta_  = v; }
    void setGamma(double v) { gamma_ = v; }

private:
    double alpha_;
    double beta_;
    double gamma_;
};

} // namespace nm