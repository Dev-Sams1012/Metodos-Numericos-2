#pragma once

#include "Matrix.hpp"
#include "core/FunctionParser.hpp"
#include <vector>
#include <string>
#include <memory>

namespace nm
{

    struct BVPResult
    {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> u;

        Matrix A;
        std::vector<double> rhs;
        int n;
        bool solved;
    };

    class BVPSolver
    {
    public:
        virtual ~BVPSolver() = default;

        BVPSolver(const BVPSolver &) = delete;
        BVPSolver &operator=(const BVPSolver &) = delete;

        virtual std::string name() const = 0;

    protected:
        BVPSolver() = default;
    };

}