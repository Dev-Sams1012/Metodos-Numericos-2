#pragma once

#include <string>
#include <stdexcept>
#include <muParser.h>

namespace nm
{

    class FunctionParser
    {
    public:
        FunctionParser();
        ~FunctionParser() = default;

        FunctionParser(const FunctionParser &) = delete;
        FunctionParser &operator=(const FunctionParser &) = delete;

        void parse(const std::string &expr);

        double eval(double x) const;

        const std::string &expression() const { return expr_; }

        bool isReady() const { return ready_; }

    private:
        mu::Parser parser_;
        mutable double xVar_{0.0};
        std::string expr_;
        bool ready_{false};
    };

}