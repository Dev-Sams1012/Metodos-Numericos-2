#ifndef derivadaLimFw_hpp
#define derivadaLimFw_hpp

#include "derivadaLimAbs.hpp"

class derivadaLimFw : public derivadaLimAbs
{
public:
    derivadaLimFw(double (*f)(double))
        : derivadaLimAbs(f, fw) {}

    double executa(double xBar, double deltaX, size_t ordem) override
    {
        if (ordem == 0)
            return f(xBar);
        else if (ordem == 1)
            return derivada1(xBar, deltaX);
        else if (ordem == 2)
            return derivada2(xBar, deltaX);

        double l = executa(xBar + deltaX, deltaX, ordem - 1);
        double r = executa(xBar, deltaX, ordem - 1);

        return (l - r) / deltaX;
    }

private:
    double derivada1(double xBar, double deltaX)
    {
        return (f(xBar + deltaX) - f(xBar)) / deltaX;
    }

    double derivada2(double xBar, double deltaX)
    {
        return (f(xBar + 2 * deltaX) - 2 * f(xBar + deltaX) + f(xBar)) / (deltaX * deltaX);
    }
};

#endif