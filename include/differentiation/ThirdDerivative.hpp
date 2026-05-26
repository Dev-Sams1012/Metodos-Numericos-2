#pragma once
// =============================================================================
// ThirdDerivative.hpp
// Aproximações por diferenças finitas para f'''(xi):
//
//   Forward O(h²): (1/h³) * [-15/2 f_i + 24 f_{i+1} - 30 f_{i+2} + 16 f_{i+3} - 3/2 f_{i+4}]
//
//   Central O(h⁴): (1 / 8h³) * [-f_{i+3} + 8f_{i+2} - 13f_{i+1} + 13f_{i-1} - 8f_{i-2} + f_{i-3}]
// =============================================================================
#include "differentiation/Differentiator.hpp"

namespace nm
{

    enum class ThirdDerivMethod
    {
        FORWARD_O2, ///< Forward  O(h²)
        CENTRAL_O4  ///< Central  O(h⁴)
    };

    class ThirdDerivative : public Differentiator
    {
    public:
        ThirdDerivative(std::shared_ptr<FunctionParser> f, ThirdDerivMethod method);

        double compute(double x, double h) const override;
        std::string name() const override;

        void setMethod(ThirdDerivMethod m) { method_ = m; }

    private:
        ThirdDerivMethod method_;

        // Métodos auxiliares privados para organização do arquivo .cpp
        double forwardO2(double x, double h) const;
        double centralO4(double x, double h) const;
    };

} // namespace nm