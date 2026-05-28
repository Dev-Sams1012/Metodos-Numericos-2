#pragma once
// =============================================================================
// BVPSolver.hpp
// Classe base abstrata para resolvedores de Problemas de Valor de Contorno
// (PVC) pelo Método das Diferenças Finitas.
//
// Um PVC genérico 1D tem a forma:
//   ED:  L[u](x) = f(x),   x ∈ [a, b]
//   CC:  u(a) = ua,  u(b) = ub
//
// Um PVC genérico 2D tem a forma:
//   ED:  ∂²u/∂x² + ∂²u/∂y² = f(x,y),  (x,y) ∈ [xa,xb]×[ya,yb]
//   CC:  u = 0 nas quatro bordas
//
// O Método das Diferenças Finitas consiste em 4 passos:
//   1. Discretizar o domínio em N partes iguais (grade de nós).
//   2. Escrever a versão discreta da ED usando diferenças centrais.
//   3. Aplicar a ED discreta em cada nó interior (gera sistema Ax = b).
//   4. Resolver o sistema linear Ax = b (via LU).
// =============================================================================
#include "Matrix.hpp"
#include "core/FunctionParser.hpp"
#include <vector>
#include <string>
#include <memory>

namespace nm {

/// @brief Resultado de um solver de PVC.
struct BVPResult {
    std::vector<double> x;      ///< Posições dos nós (1D) ou nós-x (2D)
    std::vector<double> y;      ///< Posições dos nós-y (só 2D; vazio em 1D)
    std::vector<double> u;      ///< Valores aproximados u(xi) em 1D,
                                 ///  ou u(xi,yj) ordenados por linha em 2D
    Matrix              A;      ///< Matriz do sistema (para inspeção)
    std::vector<double> rhs;    ///< Vetor independente do sistema
    int                 n;      ///< Nº de divisões (por direção)
    bool                solved; ///< Sistema foi resolvido com sucesso
};

/// @brief Interface base para resolvedores de PVC por diferenças finitas.
class BVPSolver {
public:
    virtual ~BVPSolver() = default;

    BVPSolver(const BVPSolver&)            = delete;
    BVPSolver& operator=(const BVPSolver&) = delete;

    /// @brief Descrição textual do método.
    virtual std::string name() const = 0;

protected:
    BVPSolver() = default;
};

} // namespace nm