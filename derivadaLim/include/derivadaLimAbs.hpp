#ifndef derivadaLimAbs_hpp
#define derivadaLimAbs_hpp

#include <cstddef>

class derivadaLimAbs
{
public:
    enum tipo
    {
        fw,
        bw,
        ct
    };

    double (*f)(double);
    const tipo tipoDerivada;

    virtual double executa(double xBar, double deltaX, size_t ordem) = 0;

protected:
    derivadaLimAbs(double (*f)(double), tipo t)
        : f(f), tipoDerivada(t) {}
};

#endif