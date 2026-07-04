#pragma once

#include "Matrix.hpp"
#include <vector>
#include <string>

namespace nm
{

    struct HouseholderResult
    {
        Matrix T;
        Matrix Q;
        int steps;
    };

    class HouseholderMethod
    {
    public:
        HouseholderMethod() = default;

        HouseholderResult tridiagonalize(const Matrix &A) const;

        std::string name() const { return "Transformação de Householder (tridiagonalização)"; }
    };

}