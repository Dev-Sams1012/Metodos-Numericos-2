#ifndef derivadaLimCt_hpp
#define derivadaLimCt_hpp

#include "derivadaLimAbs.hpp"

class derivadaLimCt : public derivadaLimAbs
{
public:
    derivadaLimCt(double (*f)(double));

    double executa(double xBar, double deltaX, size_t ordem) override;

private:
    double derivada1(double xBar, double deltaX);
    double derivada2(double xBar, double deltaX);
};

#endif