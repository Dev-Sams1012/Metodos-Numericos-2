#pragma once
// =============================================================================
// FDMethod2D.hpp
// Método das Diferenças Finitas para PVCs bidimensionais (Equação de Poisson).
//
// Problema resolvido (Equação de Poisson / Laplace):
//   ∂²u/∂x² + ∂²u/∂y² = f(x,y),   (x,y) ∈ [xa,xb] × [ya,yb]
//   u(x,ya) = 0   (borda inferior)
//   u(x,yb) = 0   (borda superior)
//   u(xa,y) = 0   (borda esquerda)
//   u(xb,y) = 0   (borda direita)
//
// Discretização por diferenças centrais de 2ª ordem:
//
//   ∂²u/∂x² ≈ [u_{i+1,j} - 2u_{i,j} + u_{i-1,j}] / hx²
//   ∂²u/∂y² ≈ [u_{i,j+1} - 2u_{i,j} + u_{i,j-1}] / hy²
//
// Somando, a equação discreta no nó interior (i,j) é a "máscara de 5 pontos":
//
//   (1/hx²) u_{i-1,j}  +  (1/hy²) u_{i,j-1}
// + (-2/hx² - 2/hy²) u_{i,j}
// + (1/hx²) u_{i+1,j}  +  (1/hy²) u_{i,j+1}
//   = f(xi, yj)
//
// Para grade uniforme com hx = hy = h:
//   Célula esquerda/direita : +1/h²
//   Célula inferior/superior: +1/h²
//   Célula central           : -4/h²
//
// Os nós interiores são numerados linha a linha (varrendo y crescente):
//   k = (j-1)*(nx-1) + (i-1),   i ∈ [1,nx-1], j ∈ [1,ny-1]
//
// Gera um sistema esparso (nx-1)(ny-1) × (nx-1)(ny-1) resolvido via LU.
//
// Passos do método (conforme Aula 27):
//   1. Dividir o domínio em nx×ny partes → hx = Lx/nx, hy = Ly/ny
//   2. Versão discreta da ED (diferenças centrais, máscara de 5 pontos)
//   3. Aplicar em cada nó interior (i,j) → sistema Ak = rhs
//   4. Resolver o sistema via decomposição LU
// =============================================================================
#include "bvp/BVPSolver.hpp"
#include <functional>

namespace nm {

class FDMethod2D : public BVPSolver {
public:
    FDMethod2D() = default;

    /// @brief Resolve o PVC 2D (Equação de Poisson) por diferenças finitas.
    ///
    /// @param f    Função do lado direito f(x,y). Pode ser nullptr para f ≡ 0.
    ///             A função deve aceitar duas variáveis; aqui usamos dois
    ///             FunctionParser independentes não é viável, por isso
    ///             recebemos f como uma std::function<double(double,double)>.
    /// @param xa   Extremo esquerdo em x.
    /// @param xb   Extremo direito  em x.
    /// @param ya   Extremo inferior em y.
    /// @param yb   Extremo superior em y.
    /// @param nx   Número de divisões em x.
    /// @param ny   Número de divisões em y (default = nx para grade quadrada).
    /// @param fVal Valor constante de f(x,y) quando f == nullptr.
    BVPResult solve(std::function<double(double, double)> f,
                    double xa, double xb,
                    double ya, double yb,
                    int nx, int ny = 0,
                    double fVal = 0.0) const;

    std::string name() const override {
        return "Diferenças Finitas 2D  (máscara de 5 pontos, O(h²))";
    }
};

} // namespace nm