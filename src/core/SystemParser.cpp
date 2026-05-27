// =============================================================================
// SystemParser.cpp
// =============================================================================
#include "core/SystemParser.hpp"
#include <sstream>

namespace nm
{

void SystemParser::setup(int n, const std::vector<std::string>& exprs)
{
    if (n <= 0)
        throw std::invalid_argument("SystemParser: n deve ser > 0.");
    if ((int)exprs.size() != n)
        throw std::invalid_argument("SystemParser: número de expressões != n.");

    n_     = n;
    exprs_ = exprs;
    ready_ = false;

    // Aloca variáveis — tamanho fixo, nunca redimensionar
    sVars_.assign(n, 0.0);
    tVar_ = 0.0;

    parsers_.clear();
    parsers_.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        auto p = std::make_unique<mu::Parser>();

        // Registra variáveis de estado
        for (int j = 0; j < n; ++j)
        {
            std::string var = "s" + std::to_string(j);
            p->DefineVar(var, &sVars_[j]);
        }
        
        // Registra variável de tempo
        p->DefineVar("t", &tVar_);

        // Define e testa a expressão
        try
        {
            p->SetExpr(exprs[i]);
            
            // Teste com zeros
            for (auto& v : sVars_) v = 0.0;
            tVar_ = 0.0;
            p->Eval();
        }
        catch (const mu::Parser::exception_type& e)
        {
            // mu::string_type é std::string, então não precisamos converter o e.GetMsg()
            throw std::runtime_error(
                "SystemParser expr[" + std::to_string(i) + "] \"" +
                exprs[i] + "\": " + e.GetMsg());
        }

        parsers_.push_back(std::move(p));
    }

    ready_ = true;
}

std::vector<double> SystemParser::eval(const std::vector<double>& s, double t) const
{
    if (!ready_)
        throw std::runtime_error("SystemParser::eval — parser não inicializado.");
    if ((int)s.size() != n_)
        throw std::invalid_argument("SystemParser::eval — dimensão do estado incorreta.");

    // Atualiza variáveis compartilhadas
    for (int j = 0; j < n_; ++j) sVars_[j] = s[j];
    tVar_ = t;

    std::vector<double> result(n_);
    for (int i = 0; i < n_; ++i)
    {
        try { result[i] = parsers_[i]->Eval(); }
        catch (const mu::Parser::exception_type& e)
        {
            throw std::runtime_error(
                "SystemParser::eval componente " + std::to_string(i) +
                ": " + e.GetMsg());
        }
    }
    return result;
}

} // namespace nm