#pragma once
// =============================================================================
// FunctionParser.hpp
// Encapsula o muParser para avaliação de expressões matemáticas 1D.
// =============================================================================
#include <string>
#include <stdexcept>
#include <muParser.h>

namespace nm
{

    /// @brief Wrapper do muParser para parsing e avaliação de funções matemáticas.
    ///
    /// Suporta funções de uma variável f(x).
    /// A variável disponível é 'x'.
    class FunctionParser
    {
    public:
        FunctionParser();
        ~FunctionParser() = default;

        // Impede cópias acidentais (muParser não é trivialmente copiável)
        FunctionParser(const FunctionParser &) = delete;
        FunctionParser &operator=(const FunctionParser &) = delete;

        /// @brief Define a expressão para função 1D: f(x).
        /// @throws std::runtime_error se a expressão for inválida.
        void parse(const std::string &expr);

        /// @brief Avalia f(x).
        double eval(double x) const;

        /// @brief Retorna a expressão atual.
        const std::string &expression() const { return expr_; }

        /// @brief Verifica se a expressão foi definida.
        bool isReady() const { return ready_; }

    private:
        mu::Parser parser_;
        mutable double xVar_{0.0}; ///< Variável x linkada ao parser
        std::string expr_;
        bool ready_{false};
    };

} // namespace nm