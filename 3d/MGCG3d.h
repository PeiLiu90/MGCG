#ifndef MGCG3d_H
#define MGCG3d_H

#include "MultiGrid3d.h"
#include "math.h"
#include <iostream>
#include <iomanip>
using namespace std;

class MGCG3d
{
public:
    MGCG3d(const int & ,const int &, const int &, const int &);

    void SetCoord(const double & , const double & ,const double &, const double &, const double &, const double &);
    void SetLHS(double (*)(const double &, const double &, const double &),double (*)(const double &, const double &, const double &));
    void SetRHS(double (*)(const double &, const double &, const double &));

    void Solve();
    void SolveCG();

    Vec3d<double> x;
    Vec3d<double> r;
    Vec3d<double> z;
    Vec3d<double> p;
    Vec3d<double> Ap;

    double alpha;
    double beta;
    double rz;
    double pAp;

    MG3d sys;
    int _size;
};

MGCG3d::MGCG3d(const int & q, const int & n, const int &m, const int & l):sys(q,n,m,l)
{
    _size=q;
    x.SetSize(n*pow(2,_size-1),m*pow(2,_size-1),l*pow(2,_size-1));
    r.SetSize(n*pow(2,_size-1),m*pow(2,_size-1),l*pow(2,_size-1));
    z.SetSize(n*pow(2,_size-1),m*pow(2,_size-1),l*pow(2,_size-1));
    p.SetSize(n*pow(2,_size-1),m*pow(2,_size-1),l*pow(2,_size-1));
    Ap.SetSize(n*pow(2,_size-1),m*pow(2,_size-1),l*pow(2,_size-1));
}

void MGCG3d::SetLHS(double (*g)(const double &, const double &, const double &),double (*h)(const double &, const double &, const double &))
{
    sys.SetLHS(g,h);
}

void MGCG3d::SetRHS(double (*g)(const double &, const double &,const double &))
{
    for(int i=0;i<sys.grid[_size-1].N;i++)
    {
        for(int j=0;j<sys.grid[_size-1].M;j++)
        {
            for(int k=0;k<sys.grid[_size-1].L;k++)
            {
                r(i,j,k)=g(sys.grid[_size-1]._x(i),sys.grid[_size-1]._y(j),sys.grid[_size-1]._z(k));
            }
        }
    }
}

void MGCG3d::SetCoord(const double & xmin, const double & xmax, const double & ymin, const double & ymax, const double & zmin, const double & zmax)
{
    sys.SetCoord(xmin,xmax,ymin,ymax,zmin,zmax);
}

void MGCG3d::SolveCG()
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

void MGCG3d::Solve()
{
    sys.grid[_size-1].f=r;
    sys.grid[_size-1].PreSmooth();
    sys.FineToCoarse(1);
    cout.precision(3);
    fixed(cout);
    cout<<sys.grid[_size-1].u<<endl;
//    x=0.;
//    sys.Precondition(r,z);
//    p=z;
//    double error=1.;
//    int itr=0;
//    while(error>1.E-6)
//    {
//        itr++;
//        sys.grid[_size-1].ComputeAP(p,Ap);
//        product(p,Ap,pAp);
//        product(r,z,rz);
//        alpha=rz/pAp;
//        beta=rz;
//        x+=alpha*p;
//        r-=alpha*Ap;
//        error=norm(r);
//        if(alpha<1.E-8){
//            break;
//        }
//        cout<<itr<<"  "<<error<<"  "<<alpha<<"  "<<beta<<endl;
//        sys.Precondition(r,z);
//        product(r,z,rz);
//        beta=rz/beta;
//        p=beta*p+z;
//    }
}

#endif // MGCG3d_H
