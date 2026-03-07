#include "derivadaLimBw.hpp"

derivadaLimBw::derivadaLimBw(double (*f)(double))
    : derivadaLimAbs(f, bw) {}

double derivadaLimBw::executa(double xBar, double deltaX, size_t ordem)
{
    if (ordem == 0)
        return f(xBar);
    else if (ordem == 1)
        return derivada1(xBar, deltaX);
    else if (ordem == 2)
        return derivada2(xBar, deltaX);

    double l = executa(xBar, deltaX, ordem - 1);
    double r = executa(xBar - deltaX, deltaX, ordem - 1);

    return (l - r) / deltaX;
}

double derivadaLimBw::derivada1(double xBar, double deltaX)
{
    return (f(xBar) - f(xBar - deltaX)) / deltaX;
}

double derivadaLimBw::derivada2(double xBar, double deltaX)
{
    return (f(xBar) - 2 * f(xBar - deltaX) + f(xBar - 2 * deltaX)) / (deltaX * deltaX);
}
