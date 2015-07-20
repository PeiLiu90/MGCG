#include "MGCG3d.h"
#include "Equation.h"
#include <iostream>
using namespace std;

int main()
{
    const int N=2;
    const int p=2;

    MGCG3d test(p,N,N,N);
    test.SetCoord(0,2,0,2,0,2);
    test.SetLHS(sigma,kappa);
    test.SetRHS(source);
    test.sys.grid[p-1].SetEx(exact);
    test.Solve();
    //cout<<"FINAL ERROR IS "<<norm(test.sys.grid[p-1].ex-test.x)<<endl;





//    MG3d test(p,N,N,N);
//    test.SetCoord(0,2,0,2,0,2);
//    test.SetLHS(sigma,kappa);
//    test.SetRHS(source);
//    test.grid[p-1].SetEx(exact);
//    test.Solve();
//    cout<<norm(test.grid[p-1].ex-test.grid[p-1].u)<<endl;
    return 0;
}
