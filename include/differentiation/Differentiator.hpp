#pragma once

#include "core/FunctionParser.hpp"
#include <memory>
#include <string>

namespace nm
{

    class Differentiator
    {
    public:
        explicit Differentiator(std::shared_ptr<FunctionParser> f)
            : func_(std::move(f)) {}

        virtual ~Differentiator() = default;

        Differentiator(const Differentiator &) = delete;
        Differentiator &operator=(const Differentiator &) = delete;

        virtual double compute(double x, double h) const = 0;

        virtual std::string name() const = 0;

        void setFunction(std::shared_ptr<FunctionParser> f) { func_ = std::move(f); }

    protected:
        std::shared_ptr<FunctionParser> func_;

        double f(double x, double n, double h) const
        {
            if (!func_ || !func_->isReady())
                throw std::runtime_error("Função não definida no diferenciador.");
            return func_->eval(x + n * h);
        }
    };

}