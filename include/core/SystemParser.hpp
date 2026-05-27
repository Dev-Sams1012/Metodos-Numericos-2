#pragma once
// =============================================================================
// SystemParser.hpp
// Parser para sistemas de EDOs: avalia F: ℝⁿ × ℝ → ℝⁿ
//
// Variáveis disponíveis: s0, s1, ..., s_{n-1}  (estado)  e  t  (tempo).
// Exemplo (PVI-1 escalar, n=1): "2/3 * s0"
// Exemplo (PVI-2, n=2):         ["- g - k/m * s0",  "s0"]
// =============================================================================
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <muParser.h>

namespace nm
{

class SystemParser
{
public:
    SystemParser()  = default;
    ~SystemParser() = default;

    // Não copiável (muParser armazena ponteiros para variáveis)
    SystemParser(const SystemParser&)            = delete;
    SystemParser& operator=(const SystemParser&) = delete;

    /// @brief Configura o sistema com n expressões.
    /// @param n     Dimensão do estado.
    /// @param exprs Vetor com n strings, uma por componente de F.
    /// @throws std::runtime_error se alguma expressão for inválida.
    void setup(int n, const std::vector<std::string>& exprs);

    /// @brief Avalia F(s, t).
    /// @param s Estado atual (tamanho n).
    /// @param t Instante de tempo.
    std::vector<double> eval(const std::vector<double>& s, double t) const;

    int  size()    const { return n_; }
    bool isReady() const { return ready_; }

    const std::string& expression(int i) const { return exprs_[i]; }

private:
    int  n_{0};
    bool ready_{false};
    std::vector<std::string> exprs_;

    // Armazenamento de variáveis — NUNCA redimensionar após setup()
    mutable std::vector<double> sVars_; ///< s0 .. s_{n-1}
    mutable double tVar_{0.0};

    // Um parser por componente (unique_ptr evita problemas de cópia)
    std::vector<std::unique_ptr<mu::Parser>> parsers_;
};

} // namespace nm