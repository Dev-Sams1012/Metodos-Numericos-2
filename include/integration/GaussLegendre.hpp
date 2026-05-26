#pragma once
// =============================================================================
// GaussLegendre.hpp
// Quadratura de Gauss-Legendre para ∫_a^b f(x) dx.
//
// Mudança de variável: x(s) = (a+b)/2 + (b-a)/2 * s,  s ∈ [-1,1]
// Fórmula:            I ≈ (b-a)/2 * Σ_{k=1}^{n} w_k f(x(s_k))
//
// Tabelas de nós e pesos para n = 2, 3, 4, 5.
// =============================================================================
#include "integration/Integrator.hpp"
#include <vector>

namespace nm
{

    class GaussLegendre : public Integrator
    {
    public:
        /// @param f  Função a integrar.
        /// @param n  Número de pontos (2 ≤ n ≤ 5).
        GaussLegendre(std::shared_ptr<FunctionParser> f, int n = 5);

        double compute(double a, double b) const override;
        std::string name() const override;

        void setPoints(int n);
        int points() const { return n_; }

    private:
        int n_;
        std::vector<double> nodes_;   ///< s_k ∈ [-1,1]
        std::vector<double> weights_; ///< w_k

        void loadTable(int n);
    };

} // namespace nm