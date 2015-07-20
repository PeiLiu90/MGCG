#include "MGCG3d.h"
#include "Equation.h"
#include "time.h"
#include <iostream>
using namespace std;

int main()
{
    const int N=2;
    const int p=7;

    MGCG3d test(p,N,N,N);
    test.SetCoord(0,2,0,2,0,2);
    test.SetLHS(sigma,kappa);
    test.SetRHS(source);
    test.SetEx(exact);

    clock_t start,finish;
    start = clock();
    test.Solve();
    finish = clock();
    double dur = (double)(finish - start)/CLOCKS_PER_SEC;

    cout<<"Elapse time is "<<dur<<"\nFINAL ERROR IS "<<norm(test.sys.grid[p-1].ex-test.x)<<endl;
    return 0;
}
