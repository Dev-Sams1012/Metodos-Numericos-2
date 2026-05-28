// =============================================================================
// FDMethod1D.cpp
// Método das Diferenças Finitas para PVCs unidimensionais com CCs mistas.
//
// Resolve: α·u'' + β·u' + γ·u = f(x)  em [a, b]
// com condição Dirichlet ou Neumann em cada extremidade.
//
// ── Coeficientes da equação discreta (diferenças centrais): ──────────────────
//   cL = α/h² − β/(2h)      (multiplicador de u_{i-1})
//   cC = −2α/h² + γ         (multiplicador de u_i)
//   cR = α/h² + β/(2h)      (multiplicador de u_{i+1})
//
// ── Tratamento das condições de contorno: ────────────────────────────────────
//
//   [Dirichlet esq]  u_0 = uA  (valor fixo → absorvido no rhs do nó 1)
//   [Neumann esq]    u'(a) ≈ (u_1 − u_0)/h = duA
//                    → u_0 = u_1 − h·duA
//                    → nó 0 VIRA incógnita; a equação do nó 0 é incluída no
//                      sistema com diferença progressiva O(h):
//                        u_0  incluída como 1ª linha: (cL+cC)·u_0 + cR·u_1 + ...
//                      (ver derivação no header)
//
//   [Dirichlet dir]  u_N = uB  (valor fixo → absorvido no rhs do nó N-1)
//   [Neumann dir]    u'(b) ≈ (u_N − u_{N-1})/h = duB
//                    → u_N = u_{N-1} + h·duB
//                    → nó N VIRA incógnita; última linha do sistema é
//                        cL·u_{N-1} + (cC+cR)·u_N = f(x_N) − cR·h·duB
//
// ── Mapeamento global de incógnitas: ─────────────────────────────────────────
//   Seja iFirst e iLast os índices dos nós extremos que são INCÓGNITAS:
//     D–D: iFirst=1,  iLast=N-1   (N-1 incógnitas)
//     D–N: iFirst=1,  iLast=N     (N   incógnitas)  ← questão 5
//     N–D: iFirst=0,  iLast=N-1   (N   incógnitas)
//     N–N: iFirst=0,  iLast=N     (N+1 incógnitas)
//
//   O índice de linha no sistema para o nó i é:  row = i − iFirst
// =============================================================================
#include "bvp/FDMethod1D.hpp"
#include <stdexcept>
#include <cmath>

namespace nm
{

    BVPResult FDMethod1D::solve(const std::shared_ptr<FunctionParser> &f,
                                double a, double b,
                                const BoundaryCondition &bcA,
                                const BoundaryCondition &bcB,
                                int n) const
    {
        if (n < 2)
            throw std::invalid_argument("FDMethod1D: n deve ser >= 2.");
        if (b <= a)
            throw std::invalid_argument("FDMethod1D: b deve ser > a.");
        if (std::abs(alpha_) < 1e-14)
            throw std::invalid_argument("FDMethod1D: alpha não pode ser zero.");

        // ------------------------------------------------------------------
        // Passo 1: Grade uniforme
        // ------------------------------------------------------------------
        const double h = (b - a) / static_cast<double>(n);
        const double h2 = h * h;

        // Todos os nós (incluindo fronteiras)
        BVPResult result;
        result.n = n;
        result.x.resize(n + 1);
        result.u.resize(n + 1, 0.0);
        for (int i = 0; i <= n; ++i)
            result.x[i] = a + i * h;

        // ------------------------------------------------------------------
        // Passo 2: Coeficientes da equação discreta (diferenças centrais)
        // ------------------------------------------------------------------
        const double cL = alpha_ / h2 - beta_ / (2.0 * h);
        const double cC = -2.0 * alpha_ / h2 + gamma_;
        const double cR = alpha_ / h2 + beta_ / (2.0 * h);

        // ------------------------------------------------------------------
        // Determinação do intervalo de incógnitas [iFirst, iLast]
        // ------------------------------------------------------------------
        const int iFirst = (bcA.type == BCType::DIRICHLET) ? 1 : 0;
        const int iLast = n - 1;
        const int m = iLast - iFirst + 1; // número de incógnitas

        // ------------------------------------------------------------------
        // Passo 3: Montagem do sistema m×m
        // ------------------------------------------------------------------
        result.A = Matrix(m, m, 0.0);
        result.rhs.assign(m, 0.0);

        for (int i = iFirst; i <= iLast; ++i)
        {

            const int row = i - iFirst;
            const double xi = result.x[i];
            const double fi = (f && f->isReady()) ? f->eval(xi) : 0.0;

            result.rhs[row] = fi;

            // =========================================================
            // CASO ESPECIAL: último nó com Neumann à direita
            // =========================================================
            if (i == n - 1 && bcB.type == BCType::NEUMANN)
            {

                // vizinho esquerdo
                if (i - 1 >= iFirst)
                    result.A(row, row - 1) = cL;
                else
                    result.rhs[row] -= cL * bcA.value;

                // diagonal modificada
                result.A(row, row) = cC + cR;

                // rhs modificado
                result.rhs[row] -= cR * h * bcB.value;

                continue;
            }

            // =========================================================
            // Vizinho esquerdo
            // =========================================================
            if (i - 1 < iFirst)
            {
                result.rhs[row] -= cL * bcA.value;
            }
            else
            {
                result.A(row, row - 1) = cL;
            }

            // =========================================================
            // Diagonal
            // =========================================================
            result.A(row, row) = cC;

            // =========================================================
            // Vizinho direito
            // =========================================================
            if (i + 1 > iLast)
            {

                // Dirichlet direita
                if (bcB.type == BCType::DIRICHLET)
                    result.rhs[row] -= cR * bcB.value;
            }
            else
            {
                result.A(row, row + 1) = cR;
            }
        }

        // ------------------------------------------------------------------
        // Passo 4: Resolução do sistema linear via LU
        // ------------------------------------------------------------------
        try
        {
            std::vector<double> sol = Matrix::solve(result.A, result.rhs);

            // Preenche resultado nos nós incógnita
            for (int i = iFirst; i <= iLast; ++i)
                result.u[i] = sol[i - iFirst];

            // Reconstrói nós de fronteira fixos (Dirichlet)
            if (bcA.type == BCType::DIRICHLET)
                result.u[0] = bcA.value;

            if (bcB.type == BCType::DIRICHLET)
                result.u[n] = bcB.value;

            // >>> ADICIONE AQUI <<<
            if (bcB.type == BCType::NEUMANN)
                result.u[n] = result.u[n - 1] + h * bcB.value;
                

            result.solved = true;
        }
        catch (const std::runtime_error &)
        {
            result.solved = false;
        }

        return result;
    }

} // namespace nm