#pragma once

#include "core/FunctionParser.hpp"
#include <memory>
#include <string>
#include <vector>

namespace nm
{

    enum class SpecialGaussType
    {
        HERMITE,
        LAGUERRE,
        CHEBYSHEV
    };

    class SpecialGauss
    {
    public:
        SpecialGauss(std::shared_ptr<FunctionParser> f, SpecialGaussType type, int n = 5);

        double compute() const;
        std::string name() const;

        void setFunction(std::shared_ptr<FunctionParser> f) { func_ = std::move(f); }
        void setType(SpecialGaussType t);
        void setPoints(int n);

    private:
        std::shared_ptr<FunctionParser> func_;
        SpecialGaussType type_;
        int n_;
        std::vector<double> nodes_;
        std::vector<double> weights_;

        void loadHermite(int n);
        void loadLaguerre(int n);
        void loadChebyshev(int n);
        void loadTable();
    };

}