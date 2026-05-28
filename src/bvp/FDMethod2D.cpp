// =============================================================================
// FDMethod2D.cpp
// Método das Diferenças Finitas para PVCs bidimensionais (Equação de Poisson).
//
// Resolve: ∂²u/∂x² + ∂²u/∂y² = f(x,y)
//          com u = 0 nas quatro bordas do retângulo [xa,xb]×[ya,yb].
//
// Passo 1 — Grade uniforme:
//   hx = (xb-xa)/nx,  xi = xa + i*hx,  i = 0..nx
//   hy = (yb-ya)/ny,  yj = ya + j*hy,  j = 0..ny
//
// Passo 2 — Diferenças centrais no nó (i,j):
//   ∂²u/∂x² ≈ [u_{i-1,j} - 2u_{i,j} + u_{i+1,j}] / hx²
//   ∂²u/∂y² ≈ [u_{i,j-1} - 2u_{i,j} + u_{i,j+1}] / hy²
//
//   Equação discreta (máscara de 5 pontos):
//     cW·u_{i-1,j} + cE·u_{i+1,j} + cS·u_{i,j-1} + cN·u_{i,j+1} + cC·u_{i,j}
//     = f(xi,yj)
//
//   onde:
//     cW = cE = 1/hx²
//     cS = cN = 1/hy²
//     cC = -2/hx² - 2/hy²
//
// Passo 3 — Numeração dos nós interiores (varrimento linha a linha em y):
//   k(i,j) = (j-1)*(nx-1) + (i-1),   i ∈ [1,nx-1], j ∈ [1,ny-1]
//
//   O sistema tem ordem M = (nx-1)*(ny-1).
//
// Passo 4 — Resolução via LU (Matrix::solve).
//
// Nota de convenção de indexação:
//   No resultado, u é armazenado com o índice k(i,j) acima, ou seja,
//   varrendo j de 1..ny-1 (linhas) e i de 1..nx-1 (colunas) dentro de
//   cada linha. Os nós da fronteira (u = 0) não são armazenados.
// =============================================================================
#include "bvp/FDMethod2D.hpp"
#include <stdexcept>
#include <cmath>

namespace nm {

BVPResult FDMethod2D::solve(std::function<double(double, double)> f,
                             double xa, double xb,
                             double ya, double yb,
                             int nx, int ny,
                             double fVal) const
{
    if (ny <= 0) ny = nx;   // grade quadrada por padrão
    if (nx < 2 || ny < 2)
        throw std::invalid_argument("FDMethod2D: nx e ny devem ser >= 2.");
    if (xb <= xa || yb <= ya)
        throw std::invalid_argument("FDMethod2D: domínio inválido.");

    // ------------------------------------------------------------------
    // Passo 1: Grade uniforme
    // ------------------------------------------------------------------
    const double hx  = (xb - xa) / static_cast<double>(nx);
    const double hy  = (yb - ya) / static_cast<double>(ny);
    const double hx2 = hx * hx;
    const double hy2 = hy * hy;

    const int mx = nx - 1;   // nós interiores em x
    const int my = ny - 1;   // nós interiores em y
    const int M  = mx * my;  // total de incógnitas

    // ------------------------------------------------------------------
    // Passo 2: Coeficientes da máscara de 5 pontos
    // ------------------------------------------------------------------
    const double cW = 1.0 / hx2;             // vizinho esquerdo  (i-1,j)
    const double cE = 1.0 / hx2;             // vizinho direito   (i+1,j)
    const double cS = 1.0 / hy2;             // vizinho inferior  (i,j-1)
    const double cN = 1.0 / hy2;             // vizinho superior  (i,j+1)
    const double cC = -2.0 / hx2 - 2.0 / hy2; // nó central      (i,j)

    // Função de indexação: nó interior (i,j) → índice k no vetor de incógnitas
    // i ∈ [1,mx], j ∈ [1,my]  (base-1)
    auto K = [&](int i, int j) -> int {
        return (j - 1) * mx + (i - 1);
    };

    // ------------------------------------------------------------------
    // Passo 3: Montagem do sistema M×M
    // ------------------------------------------------------------------
    BVPResult result;
    result.n = nx;
    result.A   = Matrix(M, M, 0.0);
    result.rhs.assign(M, 0.0);

    for (int j = 1; j <= my; ++j) {
        double yj = ya + j * hy;
        for (int i = 1; i <= mx; ++i) {
            double xi = xa + i * hx;
            int    k  = K(i, j);

            // Diagonal central
            result.A(k, k) = cC;

            // Vizinho esquerdo (i-1, j)
            if (i > 1)
                result.A(k, K(i - 1, j)) += cW;
            // else: fronteira x=xa → u=0, contribui 0 ao rhs

            // Vizinho direito (i+1, j)
            if (i < mx)
                result.A(k, K(i + 1, j)) += cE;
            // else: fronteira x=xb → u=0

            // Vizinho inferior (i, j-1)
            if (j > 1)
                result.A(k, K(i, j - 1)) += cS;
            // else: fronteira y=ya → u=0

            // Vizinho superior (i, j+1)
            if (j < my)
                result.A(k, K(i, j + 1)) += cN;
            // else: fronteira y=yb → u=0

            // Lado direito: f(xi, yj)
            double fij = f ? f(xi, yj) : fVal;
            result.rhs[k] = fij;

            // Nota: como todas as CC são homogêneas (u=0 na borda),
            // nenhuma parcela de contorno precisa ser subtraída do rhs.
        }
    }

    // Armazena coordenadas dos nós interiores para exibição
    result.x.resize(M);
    result.y.resize(M);
    for (int j = 1; j <= my; ++j)
        for (int i = 1; i <= mx; ++i) {
            int k = K(i, j);
            result.x[k] = xa + i * hx;
            result.y[k] = ya + j * hy;
        }

    // ------------------------------------------------------------------
    // Passo 4: Resolução do sistema linear via LU
    // ------------------------------------------------------------------
    try {
        result.u  = Matrix::solve(result.A, result.rhs);
        result.solved = true;
    }
    catch (const std::runtime_error&) {
        result.solved = false;
    }

    return result;
}

} // namespace nm