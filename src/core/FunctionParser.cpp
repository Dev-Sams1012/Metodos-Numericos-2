#include "core/FunctionParser.hpp"
#include <sstream>

namespace nm
{

    FunctionParser::FunctionParser()
    {

        parser_.DefineVar("x", &xVar_);
    }

    void FunctionParser::parse(const std::string &expr)
    {
        try
        {
            expr_ = expr;
            ready_ = false;
            parser_.SetExpr(expr);

            xVar_ = 0.0;
            parser_.Eval();
            ready_ = true;
        }
        catch (const mu::Parser::exception_type &e)
        {
            throw std::runtime_error(std::string("Erro no parser [1D]: ") + e.GetMsg());
        }
    }

    double FunctionParser::eval(double x) const
    {
        if (!ready_)
            throw std::runtime_error("Nenhuma expressão definida no parser.");
        try
        {
            xVar_ = x;
            return parser_.Eval();
        }
        catch (const mu::Parser::exception_type &e)
        {
            throw std::runtime_error(std::string("Erro na avaliação f(x): ") + e.GetMsg());
        }
    }

}