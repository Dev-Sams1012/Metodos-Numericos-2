#pragma once
// =============================================================================
// SpecialGauss.hpp
// Quadraturas especiais de Gauss (5 pontos):
//
//   Gauss-Hermite:   ∫_{-∞}^{+∞} e^{-x²} f(x) dx  ≈  Σ w_k f(x_k)
//   Gauss-Laguerre:  ∫_{0}^{+∞}  e^{-x}  f(x) dx  ≈  Σ w_k f(x_k)
//   Gauss-Chebyshev: ∫_{-1}^{1}  f(x)/√(1-x²) dx  ≈  Σ w_k f(x_k)
//                    com x_k = cos((2k-1)π/(2n)), w_k = π/n
// =============================================================================
#include "core/FunctionParser.hpp"
#include <memory>
#include <string>
#include <vector>

namespace nm {

enum class SpecialGaussType {
    HERMITE,    ///< Peso: e^{-x²}, domínio: (-∞, +∞)
    LAGUERRE,   ///< Peso: e^{-x},  domínio: [0, +∞)
    CHEBYSHEV   ///< Peso: 1/√(1-x²), domínio: [-1, 1]
};

/// @brief Quadratura de Gauss especial (5 pontos cada).
/// ATENÇÃO: f(x) deve ser apenas a parte não-peso.
/// Ex.: para Hermite, se quiser integrar g(x) = e^{-x²}·f(x), passe apenas f(x).
class SpecialGauss {
public:
    SpecialGauss(std::shared_ptr<FunctionParser> f, SpecialGaussType type, int n = 5);

    double      compute() const;
    std::string name()    const;

    void setFunction(std::shared_ptr<FunctionParser> f) { func_ = std::move(f); }
    void setType(SpecialGaussType t);
    void setPoints(int n);

private:
    std::shared_ptr<FunctionParser> func_;
    SpecialGaussType                type_;
    int                             n_;
    std::vector<double>             nodes_;
    std::vector<double>             weights_;

    void loadHermite(int n);
    void loadLaguerre(int n);
    void loadChebyshev(int n);
    void loadTable();
};

} // namespace nm