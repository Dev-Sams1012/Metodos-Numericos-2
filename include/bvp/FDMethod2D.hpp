#pragma once

#include "bvp/BVPSolver.hpp"
#include <functional>

namespace nm
{

    class FDMethod2D : public BVPSolver
    {
    public:
        FDMethod2D() = default;

        BVPResult solve(std::function<double(double, double)> f,
                        double xa, double xb,
                        double ya, double yb,
                        int nx, int ny = 0,
                        double fVal = 0.0) const;

        std::string name() const override
        {
            return "Diferenças Finitas 2D  (máscara de 5 pontos, O(h²))";
        }
    };

}