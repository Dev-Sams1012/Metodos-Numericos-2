#include "derivadaLimCt.hpp"

derivadaLimCt::derivadaLimCt(double (*f)(double))
    : derivadaLimAbs(f, ct) {}

double derivadaLimCt::executa(double xBar, double deltaX, size_t ordem)
{
    if (ordem == 0)
        return f(xBar);
    else if (ordem == 1)
        return derivada1(xBar, deltaX);
    else if (ordem == 2)
        return derivada2(xBar, deltaX);

    double l = executa(xBar + deltaX, deltaX, ordem - 1);
    double r = executa(xBar - deltaX, deltaX, ordem - 1);

    return (l - r) / (2 * deltaX);
}

double derivadaLimCt::derivada1(double xBar, double deltaX)
{
    return (f(xBar + deltaX) - f(xBar - deltaX)) / (2 * deltaX);
}

double derivadaLimCt::derivada2(double xBar, double deltaX)
{
    return (f(xBar + deltaX) - 2 * f(xBar) + f(xBar - deltaX)) / (deltaX * deltaX);
}
