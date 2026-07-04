#pragma once

#include "core/FunctionParser.hpp"
#include <memory>
#include <string>
#include <stdexcept>

namespace nm
{

    class Integrator
    {
    public:
        explicit Integrator(std::shared_ptr<FunctionParser> f)
            : func_(std::move(f)) {}

        virtual ~Integrator() = default;

        Integrator(const Integrator &) = delete;
        Integrator &operator=(const Integrator &) = delete;

        virtual double compute(double a, double b) const = 0;

        virtual std::string name() const = 0;

        void setFunction(std::shared_ptr<FunctionParser> f) { func_ = std::move(f); }

    protected:
        std::shared_ptr<FunctionParser> func_;

        double eval(double x) const
        {
            if (!func_ || !func_->isReady())
                throw std::runtime_error("Função não definida no integrador.");
            return func_->eval(x);
        }
    };

}