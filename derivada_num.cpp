double derivada1Fw(double (*f)(double), double xBar, double deltaX)
{
    return (f(xBar + deltaX) - f(xBar)) / deltaX;
}

double derivada1Bw(double (*f)(double), double xBar, double deltaX)
{
    return (f(xBar) - f(xBar - deltaX)) / deltaX;
}

double derivada1Ct(double (*f)(double), double xBar, double deltaX)
{
    return (f(xBar + deltaX) - f(xBar - deltaX) / 2 * deltaX);
}

double derivada2Fw(double (*f)(double), double xBar, double deltaX)
{
    return (f(xBar + 2 * deltaX) - 2 * f(xBar + deltaX) + f(xBar)) / (deltaX * deltaX);
}

double derivada2Bw(double (*f)(double), double xBar, double deltaX)
{
    return (f(xBar) - 2 * f(xBar - deltaX) + f(xBar - 2 * deltaX)) / (deltaX * deltaX);
}

double derivada2Ct(double (*f)(double), double xBar, double deltaX)
{
    return (f(xBar + deltaX) - 2 * f(xBar) + f(xBar - deltaX)) / (deltaX * deltaX);
}