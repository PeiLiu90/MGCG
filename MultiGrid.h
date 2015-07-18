#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "RBSGS.h"

class MG
{
public:
    MG(const int & , const int &, const int &);
    void SetCoord(const double & , const double & ,const double &, const double &);
    void SetLHS(double (*)(const double &, const double &),double (*)(const double &, const double &));
    void SetRHS(double (*)(const double &, const double &));

    void FineToCoarse(const int & k);//from k to k-1
    void CoarseToFine(const int & k);//from k to k+1

    void Solve();
    void Precondition(const Vec2d<double> & , Vec2d<double> & );
    /* 0 is the coarsest, _size-1 is the finest*/
    RBSGS * grid;
    int _size;
};

MG::MG(const int & l, const int & n, const int &m)
{
    _size=l;
    grid=new RBSGS [l];
    for(int i=0;i<_size;i++)
    {
        grid[i].SetSize(n*pow(2,i),m*pow(2,i));
    }
}

void MG::SetCoord(const double & xmin, const double & xmax, const double & ymin, const double & ymax)
{
    for(int i=0;i<_size;i++)
    {
        grid[i].SetCoord(xmin, xmax, ymin, ymax);
    }
}

void MG::SetLHS(double (*g)(const double &, const double &),double (*h)(const double &, const double &))
{
    for(int i=0;i<_size;i++)
    {
        grid[i].SetEps(g);
        grid[i].SetK2(h);
        grid[i].ComputeD();
    }
}

void MG::SetRHS(double (*g)(const double &, const double &))
{
    grid[_size-1].SetRHS(g);
}

void MG::FineToCoarse(const int & k)
{
    grid[k-1].f(0,0)= (grid[k].r(grid[k].N-1,1) + grid[k].r(0,1)*2. + grid[k].r(1,1) +
                      grid[k].r(grid[k].N-1,0)*2. + grid[k].r(0,0)*4. + grid[k].r(1,0)*2. +
                      grid[k].r(grid[k].N-1,grid[k].M-1)+ grid[k].r(0,grid[k].M-1)*2. + grid[k].r(1,grid[k].M-1))/16.;
    for(int i=1;i<grid[k-1].N;i++)
    {
        grid[k-1].f(i,0)= (grid[k].r(2*i-1,1) + grid[k].r(2*i,1)*2. + grid[k].r(2*i+1,1) +
                          grid[k].r(2*i-1,0)*2. + grid[k].r(2*i,0)*4. + grid[k].r(2*i+1,0)*2. +
                          grid[k].r(2*i-1,grid[k].M-1)+ grid[k].r(2*i,grid[k].M-1)*2. + grid[k].r(2*i+1,grid[k].M-1))/16.;
        for(int j=1;j<grid[k-1].M;j++)
        {
            grid[k-1].f(i,j)=(grid[k].r(2*i-1,2*j+1) + grid[k].r(2*i,2*j+1)*2. + grid[k].r(2*i+1,2*j+1) +
                             grid[k].r(2*i-1,2*j)*2. + grid[k].r(2*i,2*j)*4. + grid[k].r(2*i+1,2*j)*2. +
                             grid[k].r(2*i-1,2*j-1)+ grid[k].r(2*i,2*j-1)*2. + grid[k].r(2*i+1,2*j-1))/16.;
        }
    }
    for(int j=1;j<grid[k-1].M;j++)
    {
        grid[k-1].f(0,j)= (grid[k].r(grid[k].N-1,2*j+1) + grid[k].r(0,2*j+1)*2. + grid[k].r(1,2*j+1) +
                          grid[k].r(grid[k].N-1,2*j)*2. + grid[k].r(0,2*j)*4. + grid[k].r(1,2*j)*2. +
                          grid[k].r(grid[k].N-1,2*j-1)+ grid[k].r(0,2*j-1)*2. + grid[k].r(1,2*j-1))/16.;
    }
}

void MG::CoarseToFine(const int & k)
{
    for(int i=0;i<grid[k].N-1;i++)
    {
        for(int j=0;j<grid[k].M-1;j++)
        {
            grid[k+1].u(i*2,j*2)+=grid[k].u(i,j);
            grid[k+1].u(i*2,j*2+1)+=(grid[k].u(i,j)+grid[k].u(i,j+1))/2.;
            grid[k+1].u(i*2+1,j*2)+=(grid[k].u(i+1,j)+grid[k].u(i,j))/2.;
            grid[k+1].u(i*2+1,j*2+1)+=(grid[k].u(i,j)+grid[k].u(i,j+1)+grid[k].u(i+1,j)+grid[k].u(i+1,j+1))/4.;
        }
        grid[k+1].u(i*2,grid[k+1].M-2)+=grid[k].u(i,grid[k].M-1);
        grid[k+1].u(i*2,grid[k+1].M-1)+=(grid[k].u(i,grid[k].M-1)+grid[k].u(i,0))/2.;
        grid[k+1].u(i*2+1,grid[k+1].M-2)+=(grid[k].u(i+1,grid[k].M-1)+grid[k].u(i,grid[k].M-1))/2.;
        grid[k+1].u(i*2+1,grid[k+1].M-1)+=(grid[k].u(i,grid[k].M-1)+grid[k].u(i,0)+grid[k].u(i+1,grid[k].M-1)+grid[k].u(i+1,0))/4.;
    }
    for(int j=0;j<grid[k].M-1;j++)
    {
        grid[k+1].u(grid[k+1].M-2,j*2)+=grid[k].u(grid[k].M-1,j);
        grid[k+1].u(grid[k+1].M-2,j*2+1)+=(grid[k].u(grid[k].M-1,j)+grid[k].u(grid[k].M-1,j+1))/2.;
        grid[k+1].u(grid[k+1].M-1,j*2)+=(grid[k].u(0,j)+grid[k].u(grid[k].M-1,j))/2.;
        grid[k+1].u(grid[k+1].M-1,j*2+1)+=(grid[k].u(grid[k].M-1,j)+grid[k].u(grid[k].M-1,j+1)+grid[k].u(0,j)+grid[k].u(0,j+1))/4.;
    }
    grid[k+1].u(grid[k+1].M-2,grid[k+1].N-2)+=grid[k].u(grid[k].M-1,grid[k].N-1);
    grid[k+1].u(grid[k+1].M-1,grid[k+1].N-2)+=(grid[k].u(0,grid[k].N-1)+grid[k].u(grid[k].M-1,grid[k].N-1))/2.;
    grid[k+1].u(grid[k+1].M-2,grid[k+1].N-1)+=(grid[k].u(grid[k].M-1,grid[k].N-1)+grid[k].u(grid[k].M-1,0))/2.;
    grid[k+1].u(grid[k+1].M-1,grid[k+1].N-1)+=(grid[k].u(grid[k].M-1,grid[k].N-1)+grid[k].u(0,grid[k].N-1)+
                                                grid[k].u(grid[k].M-1,0)+grid[k].u(0,0))/4.;
}

void MG::Solve()
{
    int itr;
    double err;
    for(int i=_size-1;i>0;i--)
    {
        grid[i].PreSmooth();
        FineToCoarse(i);
    }
    grid[0].Solve(itr,err);
    cout<<"0th grid itr #: "<<itr<<" ,  error: "<<err<<endl;
    for(int i=0;i<_size-1;i++)
    {
        CoarseToFine(i);
        grid[i+1].Solve(itr,err);
        cout<<i+1<<"th grid itr #: "<<itr<<" ,  error: "<<err<<endl;
    }
}

void MG::Precondition(const Vec2d<double> & b, Vec2d<double> &x)
{
    grid[_size-1].f=b;
    grid[_size-1].PreSmooth();
    for(int i=_size-1;i>0;i--)
    {
        grid[i].PreSmooth();
        FineToCoarse(i);
    }
    grid[0].PreSmooth();
    for(int i=0;i<_size-1;i++)
    {
        CoarseToFine(i);
        grid[i+1].PostSmooth();
    }
    x=grid[_size-1].u;
}
#endif // MULTIGRID_H
