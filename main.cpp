#include "derivadaLimFw.hpp"
#include "derivadaLimBw.hpp"
#include "derivadaLimCt.hpp"

#include <iostream>
#include <cmath>

using namespace std;

double f(double x)
{
    return pow(x, 3) + 2 * pow(x, 2) + 3 * x + 4;
}

int main()
{
    double xBar = 2.0;
    double deltaX = 1e-4;
    size_t ordem = 1;

    derivadaLimFw fw(f);
    derivadaLimBw bw(f);
    derivadaLimCt ct(f);

    cout << "filosofia forward: " << fw.executa(xBar, deltaX, ordem) << endl;
    cout << "filosofia backward: " << bw.executa(xBar, deltaX, ordem) << endl;
    cout << "filosofia central: " << ct.executa(xBar, deltaX, ordem) << endl;

    return 0;
}