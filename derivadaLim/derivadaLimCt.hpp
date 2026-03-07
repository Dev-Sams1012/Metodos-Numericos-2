#ifndef derivadaLimCt_hpp
#define derivadaLimCt_hpp

#include "derivadaLimAbs.hpp"

class derivadaLimCt : public derivadaLimAbs
{
public:
    derivadaLimCt(double (*f)(double))
        : derivadaLimAbs(f, ct) {}

    double executa(double xBar, double deltaX, size_t ordem) override
    {
        if (ordem == 0)
            return f(xBar);
        else if (ordem == 1)
            return derivada1(xBar, deltaX);
        else if (ordem == 2)
            return derivada2(xBar, deltaX);

        double l = executa(xBar + deltaX, deltaX, ordem - 1);
        double r = executa(xBar - deltaX, deltaX, ordem - 1);

        return (l - r) / 2 * deltaX;
    }

private:
    double derivada1(double xBar, double deltaX)
    {
        return (f(xBar + deltaX) - f(xBar - deltaX) / 2 * deltaX);
    }

    double derivada2(double xBar, double deltaX)
    {
        return (f(xBar + deltaX) - 2 * f(xBar) + f(xBar - deltaX)) / (deltaX * deltaX);
    }
};

#endif