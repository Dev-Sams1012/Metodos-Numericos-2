#pragma once
// =============================================================================
// IVPSolver.hpp
// Classe base para métodos de solução de Problemas de Valor Inicial.
//
// Formato esperado do PVI (após transformação para EDO de 1ª ordem):
//   dS/dt = F(S(t), t)
//   S(t0) = S0
//
// Todos os métodos recebem SystemParser (define F), condição inicial S0,
// instante inicial t0, instante final tEnd e passo dt.
// =============================================================================
#include "core/SystemParser.hpp"
#include <vector>
#include <string>

namespace nm
{

/// @brief Resultado de um solver de PVI.
struct IVPResult
{
    std::vector<double>              t;    ///< Instantes de tempo t0, t1, ...
    std::vector<std::vector<double>> y;   ///< y[i] = estado em t[i]
    int    steps;                          ///< Total de passos realizados
    bool   completed;                      ///< Atingiu tEnd sem erro
};

/// @brief Interface base para todos os solvers de PVI.
class IVPSolver
{
public:
    virtual ~IVPSolver() = default;

    IVPSolver(const IVPSolver&)            = delete;
    IVPSolver& operator=(const IVPSolver&) = delete;

    /// @brief Resolve o PVI.
    /// @param F      Sistema de EDOs (dimensão n).
    /// @param y0     Condição inicial S(t0) (tamanho n).
    /// @param t0     Instante inicial.
    /// @param tEnd   Instante final.
    /// @param dt     Passo de tempo.
    /// @param maxSteps Limite de passos (segurança).
    virtual IVPResult solve(const SystemParser& F,
                             const std::vector<double>& y0,
                             double t0, double tEnd, double dt,
                             int maxSteps = 200000) const = 0;

    virtual std::string name() const = 0;

protected:
    IVPSolver() = default;

    // ---- Helpers vetoriais usados pelas subclasses -------------------------

    /// @brief  y_out = a + alpha * b
    static std::vector<double> axpy(const std::vector<double>& a,
                                     double alpha,
                                     const std::vector<double>& b)
    {
        int n = (int)a.size();
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i) r[i] = a[i] + alpha * b[i];
        return r;
    }

    /// @brief Norma Euclidiana de um vetor.
    static double norm(const std::vector<double>& v)
    {
        double s = 0.0;
        for (double x : v) s += x * x;
        return std::sqrt(s);
    }

    /// @brief Diferença normalizada ||a - b|| / (1 + ||b||).
    static double relDiff(const std::vector<double>& a, const std::vector<double>& b)
    {
        int n = (int)a.size();
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (a[i] - b[i]) * (a[i] - b[i]);
        return std::sqrt(diff) / (1.0 + norm(b));
    }
};

} // namespace nm