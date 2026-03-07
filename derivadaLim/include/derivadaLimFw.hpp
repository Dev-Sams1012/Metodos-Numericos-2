#ifndef derivadaLimFw_hpp
#define derivadaLimFw_hpp

#include "derivadaLimAbs.hpp"

class derivadaLimFw : public derivadaLimAbs
{
public:
    derivadaLimFw(double (*f)(double));

    double executa(double xBar, double deltaX, size_t ordem) override;

private:
    double derivada1(double xBar, double deltaX);
    double derivada2(double xBar, double deltaX);
};

#endif