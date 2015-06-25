#include "Helmholtz2d.h"

int main()
{
    const int Nx=3;
    const int Ny=3;
    const double h = 2./double(Nx);
    const int M=2;
    Helmholtz2d<double,Nx,Ny,M> test(h);
    test.SetMat();
    test.SetRHS();
    test.Solve();
    return 0;
}
