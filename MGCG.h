#ifndef MGCG_H
#define MGCG_H

#include "MultiGrid.h"
#include "math.h"
#include <iostream>
#include <iomanip>
using namespace std;

class MGCG
{
public:
    MGCG(const int & ,const int &, const int &);

    void SetCoord(const double & , const double & ,const double &, const double &);
    void SetLHS(double (*)(const double &, const double &),double (*)(const double &, const double &));
    void SetRHS(double (*)(const double &, const double &));

    void Solve();
    void SolveCG();

    Vec2d<double> x;
    Vec2d<double> r;
    Vec2d<double> z;
    Vec2d<double> p;
    Vec2d<double> Ap;

    double alpha;
    double beta;
    double rz;
    double pAp;

    MG sys;
    int _size;
};

MGCG::MGCG(const int & l, const int & n, const int &m):sys(l,n,m)
{
    _size=l;
    x.SetSize(n*pow(2,_size-1),m*pow(2,_size-1));
    r.SetSize(n*pow(2,_size-1),m*pow(2,_size-1));
    z.SetSize(n*pow(2,_size-1),m*pow(2,_size-1));
    p.SetSize(n*pow(2,_size-1),m*pow(2,_size-1));
    Ap.SetSize(n*pow(2,_size-1),m*pow(2,_size-1));
}

void MGCG::SetLHS(double (*g)(const double &, const double &),double (*h)(const double &, const double &))
{
    sys.SetLHS(g,h);
}

void MGCG::SetRHS(double (*g)(const double &, const double &))
{
    for(int i=0;i<sys.grid[_size-1].N;i++)
    {
        for(int j=0;j<sys.grid[_size-1].M;j++)
        {
            r(i,j)=g(sys.grid[_size-1]._x(i),sys.grid[_size-1]._y(j));
        }
    }
}

void MGCG::SetCoord(const double & xmin, const double & xmax, const double & ymin, const double & ymax)
{
    sys.SetCoord(xmin,xmax,ymin,ymax);
}

void MGCG::SolveCG()
{
    x=0.;
    p=r;
    double error=1.;
    product(r,r,rz);
    int itr=0;
    while(error>1.E-6)
    {
        itr++;
        sys.grid[_size-1].ComputeAP(p,Ap);
        product(p,Ap,pAp);
        alpha=rz/pAp;
        x+=alpha*p;
        r-=alpha*Ap;
        error=norm(r);
        cout<<itr<<"  "<<error<<endl;
        beta=rz;
        product(r,r,rz);
        beta=rz/beta;
        p=beta*p+r;
    }
}

void MGCG::Solve()
{
    x=0.;
    sys.Precondition(r,z);
    p=z;
    double error=1.;
    int itr=0;
    while(error>1.E-6)
    {
        itr++;
        sys.grid[_size-1].ComputeAP(p,Ap);
        product(p,Ap,pAp);
        product(r,z,rz);
        alpha=rz/pAp;
        beta=rz;
        x+=alpha*p;
        r-=alpha*Ap;
        error=norm(r);
        if(alpha<1.E-8){
            break;
        }
        cout<<itr<<"  "<<error<<"  "<<alpha<<"  "<<beta<<endl;
        sys.Precondition(r,z);
        product(r,z,rz);
        beta=rz/beta;
        p=beta*p+z;
    }
}

#endif // MGCG_H
