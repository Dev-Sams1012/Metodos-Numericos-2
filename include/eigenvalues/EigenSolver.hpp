#pragma once
// =============================================================================
// EigenSolver.hpp
// Classe base abstrata para resolvedores de autovalores/autovetores.
// Todos os métodos operam sobre matrizes simétricas reais.
// =============================================================================
#include "Matrix.hpp"
#include <vector>
#include <string>
#include <stdexcept>

namespace nm {

/// @brief Resultado de um método de autovalores.
struct EigenResult {
    double              eigenvalue;   ///< Autovalor dominante (ou alvo)
    std::vector<double> eigenvector;  ///< Autovetor correspondente
    int                 iterations;   ///< Iterações até convergência
    bool                converged;    ///< Critério de convergência satisfeito?
};

/// @brief Interface base para solvers de autovalores.
class EigenSolver {
public:
    virtual ~EigenSolver() = default;

    EigenSolver(const EigenSolver&)            = delete;
    EigenSolver& operator=(const EigenSolver&) = delete;

    /// @brief Executa o algoritmo.
    /// @param A       Matriz simétrica de entrada.
    /// @param tol     Tolerância de convergência (ε).
    /// @param maxIter Número máximo de iterações.
    virtual EigenResult solve(const Matrix& A, double tol = 1e-10,
                              int maxIter = 1000) const = 0;

    /// @brief Descrição do método.
    virtual std::string name() const = 0;

protected:
    EigenSolver() = default;
};

} // namespace nm