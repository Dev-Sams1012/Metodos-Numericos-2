#pragma once
// =============================================================================
// Differentiator.hpp
// Classe base abstrata para todos os métodos de diferenciação numérica.
// Segue o padrão Template Method / Strategy — subclasses implementam compute().
// =============================================================================
#include "core/FunctionParser.hpp"
#include <memory>
#include <string>

namespace nm
{

    /// @brief Interface base para diferenciadores numéricos.
    class Differentiator
    {
    public:
        explicit Differentiator(std::shared_ptr<FunctionParser> f)
            : func_(std::move(f)) {}

        virtual ~Differentiator() = default;

        // Não copiável para evitar aliasing do FunctionParser
        Differentiator(const Differentiator &) = delete;
        Differentiator &operator=(const Differentiator &) = delete;

        /// @brief Calcula a derivada em x usando passo h.
        /// @param x  Ponto de avaliação.
        /// @param h  Tamanho do passo (deve ser > 0).
        /// @returns  Aproximação da derivada.
        virtual double compute(double x, double h) const = 0;

        /// @brief Retorna descrição textual do método.
        virtual std::string name() const = 0;

        /// @brief Substitui a função a ser diferenciada.
        void setFunction(std::shared_ptr<FunctionParser> f) { func_ = std::move(f); }

    protected:
        std::shared_ptr<FunctionParser> func_;

        /// @brief Avalia f(xi + n*h) com verificação de ponteiro.
        double f(double x, double n, double h) const
        {
            if (!func_ || !func_->isReady())
                throw std::runtime_error("Função não definida no diferenciador.");
            return func_->eval(x + n * h);
        }
    };

} // namespace nm