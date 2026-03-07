#ifndef derivadaLimBw_hpp
#define derivadaLimBw_hpp

#include "derivadaLimAbs.hpp"

class derivadaLimBw : public derivadaLimAbs
{
public:
    derivadaLimBw(double (*f)(double));

    double executa(double xBar, double deltaX, size_t ordem) override;

private:
    double derivada1(double xBar, double deltaX);
    double derivada2(double xBar, double deltaX);
};

#endif