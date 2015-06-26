#include "Helmholtz2d.h"

int main()
{
    const int Nx=4;
    const int Ny=4;
    const double h = 2./double(Nx);
    const int M=7;
    Helmholtz2d<double,Nx,Ny,M> test(h);
    test.SetMat();
    test.SetRHS();
    test.Solve();

    cout<<norm(test._unknown[M-1]-test._exact)<<endl;

    Vec2d<double> rhs(Nx*64,Ny*64);
    test._mat[M-1].Product(test._unknown[M-1],rhs);
    Vec2d<double> rhs2(Nx*64,Ny*64);
    test._mat[M-1].Product(test._exact,rhs2);
    cout<<norm(rhs-rhs2)<<"  "<<norm(rhs-test._b)<<"  "<<norm(rhs2-test._b)<<endl;
    return 0;
}
