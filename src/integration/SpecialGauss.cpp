// =============================================================================
// SpecialGauss.cpp
// Tabelas de Gauss-Hermite e Gauss-Laguerre (5 pontos).
// Referências: Abramowitz & Stegun, cap. 25.
// =============================================================================
#include "integration/SpecialGauss.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace nm {

SpecialGauss::SpecialGauss(std::shared_ptr<FunctionParser> f,
                           SpecialGaussType type, int n)
    : func_(std::move(f)), type_(type), n_(n) {
    loadTable();
}

void SpecialGauss::setType(SpecialGaussType t) { type_ = t; loadTable(); }
void SpecialGauss::setPoints(int n)            { n_ = n;   loadTable(); }

void SpecialGauss::loadTable() {
    switch (type_) {
        case SpecialGaussType::HERMITE:   loadHermite(n_);   break;
        case SpecialGaussType::LAGUERRE:  loadLaguerre(n_);  break;
        case SpecialGaussType::CHEBYSHEV: loadChebyshev(n_); break;
    }
}

// ---------- Gauss-Hermite -------------------------------------------------
// Tabela de 5 pontos para ∫_{-∞}^{∞} e^{-x²} f(x) dx
void SpecialGauss::loadHermite(int n) {
    if (n != 5)
        throw std::invalid_argument("Gauss-Hermite: apenas n=5 implementado.");
    // Raízes do polinômio de Hermite H_5(x)
    nodes_ = {
        -2.0201828704560856, -0.9585724646137554,
         0.0,
         0.9585724646137554,  2.0201828704560856
    };
    // Pesos correspondentes
    weights_ = {
        0.019953242059045923, 0.39361932315224975,
        0.9453087204829419,
        0.39361932315224975,  0.019953242059045923
    };
}

// ---------- Gauss-Laguerre ------------------------------------------------
// Tabela de 5 pontos para ∫_0^{∞} e^{-x} f(x) dx
void SpecialGauss::loadLaguerre(int n) {
    if (n != 5)
        throw std::invalid_argument("Gauss-Laguerre: apenas n=5 implementado.");
    // Raízes do polinômio de Laguerre L_5(x)
    nodes_ = {
        0.26356031971814109, 1.4134030591065167,
        3.5964257710407243,  7.0858100058588376,
        12.640800844275783
    };
    // Pesos
    weights_ = {
        0.52175561058280873, 0.39866681108315199,
        0.075942449681707595, 0.0036117586799220484,
        0.00002336997238577622
    };
}

// ---------- Gauss-Chebyshev do 1o. tipo ----------------------------------
// ∫_{-1}^{1} f(x)/√(1-x²) dx  ≈  (π/n) Σ f(x_k)
// x_k = cos((2k-1)π/(2n))
void SpecialGauss::loadChebyshev(int n) {
    if (n < 1) throw std::invalid_argument("Gauss-Chebyshev: n deve ser >= 1.");
    nodes_.resize(n);
    weights_.resize(n);
    const double pi = std::acos(-1.0);
    for (int k = 1; k <= n; ++k) {
        nodes_[k - 1]   = std::cos((2.0 * k - 1.0) * pi / (2.0 * n));
        weights_[k - 1] = pi / n;
    }
}

double SpecialGauss::compute() const {
    if (!func_ || !func_->isReady())
        throw std::runtime_error("SpecialGauss: função não definida.");
    double result = 0.0;
    for (int k = 0; k < (int)nodes_.size(); ++k)
        result += weights_[k] * func_->eval(nodes_[k]);
    return result;
}

std::string SpecialGauss::name() const {
    std::ostringstream oss;
    switch (type_) {
        case SpecialGaussType::HERMITE:
            oss << "Gauss-Hermite   (n=" << n_ << ", peso e^{-x²})";   break;
        case SpecialGaussType::LAGUERRE:
            oss << "Gauss-Laguerre  (n=" << n_ << ", peso e^{-x})";    break;
        case SpecialGaussType::CHEBYSHEV:
            oss << "Gauss-Chebyshev (n=" << n_ << ", peso 1/√(1-x²))"; break;
    }
    return oss.str();
}

} // namespace nm