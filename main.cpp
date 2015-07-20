#include "MGCG.h"
#include "Equation.h"
#include <iostream>
#include <iomanip>
using namespace std;

int main()
{
    const int L=9;
    const int N=2;
    const int M=N;
    MGCG test(L,N,M);
    test.SetCoord(0,2,0,2);
    test.SetLHS(sigma,kappa);
    test.SetRHS(source);
    test.Solve();
    return 0;
}





//int main()
//{
//    const int N=100;
//    const int M=N;
//
//    RBSGS test;
//    test.SetSize(N,M);
//    test.SetCoord(0,2,0,2);
//    test.SetEps(sigma);
//    test.SetK2(kappa);
//    test.ComputeD();
//    test.SetRHS(source);
//    Vec2d<double> ex(N,M);
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<M;j++)
//        {
//            ex(i,j)=exact(test._x(i),test._y(j));
//        }
//    }
//    test.Solve();
//    cout<<norm(ex-test.u);
//    return 0;
//}
