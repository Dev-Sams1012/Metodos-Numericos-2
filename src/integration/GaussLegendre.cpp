// =============================================================================
// GaussLegendre.cpp
// Tabelas de nós e pesos de Gauss-Legendre para n = 2 a 5.
// Referência: Abramowitz & Stegun, Tabela 25.4.
// =============================================================================
#include "integration/GaussLegendre.hpp"
#include <stdexcept>
#include <sstream>

namespace nm
{

    GaussLegendre::GaussLegendre(std::shared_ptr<FunctionParser> f, int n)
        : Integrator(std::move(f))
    {
        setPoints(n);
    }

    void GaussLegendre::setPoints(int n)
    {
        if (n < 2 || n > 5)
            throw std::invalid_argument("GaussLegendre: n deve estar entre 2 e 5.");
        n_ = n;
        loadTable(n);
    }

    void GaussLegendre::loadTable(int n)
    {
        nodes_.clear();
        weights_.clear();

        switch (n)
        {
        case 2:
            // n=2: nós ±1/√3, pesos 1,1
            nodes_ = {-0.5773502691896257, 0.5773502691896257};
            weights_ = {1.0, 1.0};
            break;

        case 3:
            // n=3: nós 0, ±√(3/5); pesos 8/9, 5/9
            nodes_ = {-0.7745966692414834, 0.0, 0.7745966692414834};
            weights_ = {0.5555555555555556, 0.8888888888888888,
                        0.5555555555555556};
            break;

        case 4:
            // n=4: nós ±0.3399810..., ±0.8611363...
            // pesos 0.6521451..., 0.3478548...
            nodes_ = {
                -0.8611363115940526, -0.3399810435848563,
                0.3399810435848563, 0.8611363115940526};
            weights_ = {
                0.3478548451374538, 0.6521451548625461,
                0.6521451548625461, 0.3478548451374538};
            break;

        case 5:
            // n=5: nós 0, ±0.5384693..., ±0.9061798...
            // pesos 0.5688888..., 0.4786286..., 0.2369268...
            nodes_ = {
                -0.9061798459386640, -0.5384693101056831,
                0.0,
                0.5384693101056831, 0.9061798459386640};
            weights_ = {
                0.2369268850561891, 0.4786286704993665,
                0.5688888888888889,
                0.4786286704993665, 0.2369268850561891};
            break;
        }
    }

    double GaussLegendre::compute(double a, double b) const
    {
        // Mudança de variável: x(s) = (a+b)/2 + (b-a)/2 * s
        // I ≈ (b-a)/2 * Σ w_k * f(x(s_k))
        const double mid = (a + b) / 2.0;
        const double half = (b - a) / 2.0;
        double result = 0.0;
        for (int k = 0; k < n_; ++k)
        {
            double xk = mid + half * nodes_[k];
            result += weights_[k] * eval(xk);
        }
        return half * result;
    }

    std::string GaussLegendre::name() const
    {
        std::ostringstream oss;
        oss << "Gauss-Legendre (n=" << n_ << ")";
        return oss.str();
    }

} // namespace nm