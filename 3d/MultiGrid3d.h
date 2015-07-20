#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "RBSGS3d.h"

class MG3d
{
public:
    MG3d(const int & , const int &, const int &, const int & );
    void SetCoord(const double & , const double & ,const double &, const double &, const double &, const double &);
    void SetLHS(double (*)(const double &, const double &, const double &),double (*)(const double &, const double &, const double &));
    void SetRHS(double (*)(const double &, const double &, const double &));

    void FineToCoarse(const int & k);//from k to k-1
    void CoarseToFine(const int & k);//from k to k+1

    void Solve();
    void Precondition(const Vec3d<double> & , Vec3d<double> & );
    /* 0 is the coarsest, _size-1 is the finest*/
    RBSGS3d * grid;
    int _size;
};

MG3d::MG3d(const int & p, const int & n, const int &m, const int & l)
{
    _size=p;
    grid=new RBSGS3d [p];
    for(int i=0;i<_size;i++)
    {
        grid[i].SetSize(n*pow(2,i),m*pow(2,i),l*pow(2,i));
    }
}

void MG3d::SetCoord(const double & xmin, const double & xmax, const double & ymin, const double & ymax, const double & zmin, const double & zmax)
{
    for(int i=0;i<_size;i++)
    {
        grid[i].SetCoord(xmin, xmax, ymin, ymax, zmin, zmax);
    }
}

void MG3d::SetLHS(double (*g)(const double &, const double &, const double &),double (*h)(const double &, const double &, const double &))
{
    for(int i=0;i<_size;i++)
    {
        grid[i].SetEps(g);
        grid[i].SetK2(h);
        grid[i].ComputeD();
    }
}

void MG3d::SetRHS(double (*g)(const double &, const double &, const double &))
{
    grid[_size-1].SetRHS(g);
}

void MG3d::FineToCoarse(const int & p)
{
    grid[p-1].f(0,0,0) = (grid[p].r(0,0,0)*8.
                             +grid[p].r(1,0,0)*4. + grid[p].r(grid[p].N-1,0,0)*4.
                             +grid[p].r(0,1,0)*4. + grid[p].r(0,grid[p].M-1,0)*4.
                             +grid[p].r(0,0,1)*4. + grid[p].r(0,0,grid[p].L-1)*4.
                             +grid[p].r(1,1,0)*2. +grid[p].r(1,grid[p].M-1,0)*2. + grid[p].r(grid[p].N-1,1,0)*2.+ grid[p].r(grid[p].N-1,grid[p].M-1,0)*2.
                             +grid[p].r(1,0,1)*2. +grid[p].r(1,0,grid[p].L-1)*2. + grid[p].r(grid[p].N-1,0,1)*2.+ grid[p].r(grid[p].N-1,0,grid[p].L-1)*2.
                             +grid[p].r(0,1,1)*2. +grid[p].r(0,1,grid[p].L-1)*2. + grid[p].r(0,grid[p].M-1,1)*2.+ grid[p].r(0,grid[p].M-1,grid[p].L-1)*2.
                             +grid[p].r(1,1,1) +grid[p].r(grid[p].N-1,1,1)+grid[p].r(1,grid[p].M-1,1)+grid[p].r(1,1,grid[p].L-1)+grid[p].r(grid[p].N-1,grid[p].M-1,1)+grid[p].r(grid[p].N-1,1,grid[p].L-1)+grid[p].r(1,grid[p].M-1,grid[p].L-1)+grid[p].r(grid[p].N-1,grid[p].M-1,grid[p].L-1))/64.;
    for(int i=1;i<grid[p-1].N;i++)
    {
        grid[p-1].f(i,0,0) = (grid[p].r(2*i,0,0)*8.
                             +grid[p].r(2*i+1,0,0)*4. + grid[p].r(2*i-1,0,0)*4.
                             +grid[p].r(2*i,1,0)*4. + grid[p].r(2*i,grid[p].M-1,0)*4.
                             +grid[p].r(2*i,0,1)*4. + grid[p].r(2*i,0,grid[p].L-1)*4.
                             +grid[p].r(2*i+1,1,0)*2. +grid[p].r(2*i+1,grid[p].M-1,0)*2. + grid[p].r(2*i-1,1,0)*2.+ grid[p].r(2*i-1,grid[p].M-1,0)*2.
                             +grid[p].r(2*i+1,0,1)*2. +grid[p].r(2*i+1,0,grid[p].L-1)*2. + grid[p].r(2*i-1,0,1)*2.+ grid[p].r(2*i-1,0,grid[p].L-1)*2.
                             +grid[p].r(2*i,1,1)*2. +grid[p].r(2*i,1,grid[p].L-1)*2. + grid[p].r(2*i,grid[p].M-1,1)*2.+ grid[p].r(2*i,grid[p].M-1,grid[p].L-1)*2.
                             +grid[p].r(2*i+1,1,1) +grid[p].r(2*i-1,1,1)+grid[p].r(2*i+1,grid[p].M-1,1)+grid[p].r(2*i+1,1,grid[p].L-1)+grid[p].r(2*i-1,grid[p].M-1,1)+grid[p].r(2*i-1,1,grid[p].L-1)+grid[p].r(2*i+1,grid[p].M-1,grid[p].L-1)+grid[p].r(2*i-1,grid[p].M-1,grid[p].L-1))/64.;
        for(int j=1;j<grid[p-1].M;j++)
        {
            grid[p-1].f(i,j,0) = (grid[p].r(2*i,2*j,0)*8.
                                 +grid[p].r(2*i+1,2*j,0)*4. + grid[p].r(2*i-1,2*j,0)*4.
                                 +grid[p].r(2*i,2*j+1,0)*4. + grid[p].r(2*i,2*j-1,0)*4.
                                 +grid[p].r(2*i,2*j,1)*4. + grid[p].r(2*i,2*j,grid[p].L-1)*4.
                                 +grid[p].r(2*i+1,2*j+1,0)*2. +grid[p].r(2*i+1,2*j-1,0)*2. + grid[p].r(2*i-1,2*j+1,0)*2.+ grid[p].r(2*i-1,2*j-1,0)*2.
                                 +grid[p].r(2*i+1,2*j,1)*2. +grid[p].r(2*i+1,2*j,grid[p].L-1)*2. + grid[p].r(2*i-1,2*j,1)*2.+ grid[p].r(2*i-1,2*j,grid[p].L-1)*2.
                                 +grid[p].r(2*i,2*j+1,1)*2. +grid[p].r(2*i,2*j+1,grid[p].L-1)*2. + grid[p].r(2*i,2*j-1,1)*2.+ grid[p].r(2*i,2*j-1,grid[p].L-1)*2.
                                 +grid[p].r(2*i+1,2*j+1,1) +grid[p].r(2*i-1,2*j+1,1)+grid[p].r(2*i+1,2*j-1,1)+grid[p].r(2*i+1,2*j+1,grid[p].L-1)+grid[p].r(2*i-1,2*j-1,1)+grid[p].r(2*i-1,2*j+1,grid[p].L-1)+grid[p].r(2*i+1,2*j-1,grid[p].L-1)+grid[p].r(2*i-1,2*j-1,grid[p].L-1))/64.;
            for(int k=1;k<grid[p-1].L;k++)
            {
                grid[p-1].f(i,j,k) = (grid[p].r(2*i,2*j,2*k)*8.
                                     +grid[p].r(2*i+1,2*j,2*k)*4. + grid[p].r(2*i-1,2*j,2*k)*4.
                                     +grid[p].r(2*i,2*j+1,2*k)*4. + grid[p].r(2*i,2*j-1,2*k)*4.
                                     +grid[p].r(2*i,2*j,2*k+1)*4. + grid[p].r(2*i,2*j,2*k-1)*4.
                                     +grid[p].r(2*i+1,2*j+1,2*k)*2. +grid[p].r(2*i+1,2*j-1,2*k)*2. + grid[p].r(2*i-1,2*j+1,2*k)*2.+ grid[p].r(2*i-1,2*j-1,2*k)*2.
                                     +grid[p].r(2*i+1,2*j,2*k+1)*2. +grid[p].r(2*i+1,2*j,2*k-1)*2. + grid[p].r(2*i-1,2*j,2*k+1)*2.+ grid[p].r(2*i-1,2*j,2*k-1)*2.
                                     +grid[p].r(2*i,2*j+1,2*k+1)*2. +grid[p].r(2*i,2*j+1,2*k-1)*2. + grid[p].r(2*i,2*j-1,2*k+1)*2.+ grid[p].r(2*i,2*j-1,2*k-1)*2.
                                     +grid[p].r(2*i+1,2*j+1,2*k+1) +grid[p].r(2*i-1,2*j+1,2*k+1)+grid[p].r(2*i+1,2*j-1,2*k+1)+grid[p].r(2*i+1,2*j+1,2*k-1)+grid[p].r(2*i-1,2*j-1,2*k+1)+grid[p].r(2*i-1,2*j+1,2*k-1)+grid[p].r(2*i+1,2*j-1,2*k-1)+grid[p].r(2*i-1,2*j-1,2*k-1))/64.;
            }
        }
        for(int k=1;k<grid[p-1].L;k++)
        {
            grid[p-1].f(i,0,k) =(grid[p].r(2*i,0,2*k)*8.
                             +grid[p].r(2*i+1,0,2*k)*4. + grid[p].r(2*i-1,0,2*k)*4.
                             +grid[p].r(2*i,1,2*k)*4. + grid[p].r(2*i,grid[p].M-1,2*k)*4.
                             +grid[p].r(2*i,0,2*k+1)*4. + grid[p].r(2*i,0,2*k-1)*4.
                             +grid[p].r(2*i+1,1,2*k)*2. +grid[p].r(2*i+1,grid[p].M-1,2*k)*2. + grid[p].r(2*i-1,1,2*k)*2.+ grid[p].r(2*i-1,grid[p].M-1,2*k)*2.
                             +grid[p].r(2*i+1,0,2*k+1)*2. +grid[p].r(2*i+1,0,2*k-1)*2. + grid[p].r(2*i-1,0,2*k+1)*2.+ grid[p].r(2*i-1,0,2*k-1)*2.
                             +grid[p].r(2*i,1,2*k+1)*2. +grid[p].r(2*i,1,2*k-1)*2. + grid[p].r(2*i,grid[p].M-1,2*k+1)*2.+ grid[p].r(2*i,grid[p].M-1,2*k-1)*2.
                             +grid[p].r(2*i+1,1,2*k+1) +grid[p].r(2*i-1,1,2*k+1)+grid[p].r(2*i+1,grid[p].M-1,2*k+1)+grid[p].r(2*i+1,1,2*k-1)+grid[p].r(2*i-1,grid[p].M-1,2*k+1)+grid[p].r(2*i-1,1,2*k-1)+grid[p].r(2*i+1,grid[p].M-1,2*k-1)+grid[p].r(2*i-1,grid[p].M-1,2*k-1))/64.;
        }
    }
    for(int j=1;j<grid[p-1].M;j++)
    {
        grid[p-1].f(0,j,0)=(grid[p].r(0,2*j,0)*8.
                         +grid[p].r(1,2*j,0)*4. + grid[p].r(grid[p].N-1,2*j,0)*4.
                         +grid[p].r(0,2*j+1,0)*4. + grid[p].r(0,2*j-1,0)*4.
                         +grid[p].r(0,2*j,1)*4. + grid[p].r(0,2*j,grid[p].L-1)*4.
                         +grid[p].r(1,2*j+1,0)*2. +grid[p].r(1,2*j-1,0)*2. + grid[p].r(grid[p].N-1,2*j+1,0)*2.+ grid[p].r(grid[p].N-1,2*j-1,0)*2.
                         +grid[p].r(1,2*j,1)*2. +grid[p].r(1,2*j,grid[p].L-1)*2. + grid[p].r(grid[p].N-1,2*j,1)*2.+ grid[p].r(grid[p].N-1,2*j,grid[p].L-1)*2.
                         +grid[p].r(0,2*j+1,1)*2. +grid[p].r(0,2*j+1,grid[p].L-1)*2. + grid[p].r(0,2*j-1,1)*2.+ grid[p].r(0,2*j-1,grid[p].L-1)*2.
                         +grid[p].r(1,2*j+1,1) +grid[p].r(grid[p].N-1,2*j+1,1)+grid[p].r(1,2*j-1,1)+grid[p].r(1,2*j+1,grid[p].L-1)+grid[p].r(grid[p].N-1,2*j-1,1)+grid[p].r(grid[p].N-1,2*j+1,grid[p].L-1)+grid[p].r(1,2*j-1,grid[p].L-1)+grid[p].r(grid[p].N-1,2*j-1,grid[p].L-1))/64.;
        for(int k=1;k<grid[p-1].L;k++)
        {
            grid[p-1].f(0,j,k)=(grid[p].r(0,2*j,2*k)*8.
                             +grid[p].r(1,2*j,2*k)*4. + grid[p].r(grid[p].N-1,2*j,2*k)*4.
                             +grid[p].r(0,2*j+1,2*k)*4. + grid[p].r(0,2*j-1,2*k)*4.
                             +grid[p].r(0,2*j,2*k+1)*4. + grid[p].r(0,2*j,2*k-1)*4.
                             +grid[p].r(1,2*j+1,2*k)*2. +grid[p].r(1,2*j-1,2*k)*2. + grid[p].r(grid[p].N-1,2*j+1,2*k)*2.+ grid[p].r(grid[p].N-1,2*j-1,2*k)*2.
                             +grid[p].r(1,2*j,2*k+1)*2. +grid[p].r(1,2*j,2*k-1)*2. + grid[p].r(grid[p].N-1,2*j,2*k+1)*2.+ grid[p].r(grid[p].N-1,2*j,2*k-1)*2.
                             +grid[p].r(0,2*j+1,2*k+1)*2. +grid[p].r(0,2*j+1,2*k-1)*2. + grid[p].r(0,2*j-1,2*k+1)*2.+ grid[p].r(0,2*j-1,2*k-1)*2.
                             +grid[p].r(1,2*j+1,2*k+1) +grid[p].r(grid[p].N-1,2*j+1,2*k+1)+grid[p].r(1,2*j-1,2*k+1)+grid[p].r(1,2*j+1,2*k-1)+grid[p].r(grid[p].N-1,2*j-1,2*k+1)+grid[p].r(grid[p].N-1,2*j+1,2*k-1)+grid[p].r(1,2*j-1,2*k-1)+grid[p].r(grid[p].N-1,2*j-1,2*k-1))/64.;
        }
    }
    for(int k=1;k<grid[p-1].L;k++)
    {
        grid[p-1].f(0,0,k)=(grid[p].r(0,0,2*k)*8.
                         +grid[p].r(1,0,2*k)*4. + grid[p].r(grid[p].N-1,0,2*k)*4.
                         +grid[p].r(0,1,2*k)*4. + grid[p].r(0,grid[p].M-1,2*k)*4.
                         +grid[p].r(0,0,2*k+1)*4. + grid[p].r(0,0,2*k-1)*4.
                         +grid[p].r(1,1,2*k)*2. +grid[p].r(1,grid[p].M-1,2*k)*2. + grid[p].r(grid[p].N-1,1,2*k)*2.+ grid[p].r(grid[p].N-1,grid[p].M-1,2*k)*2.
                         +grid[p].r(1,0,2*k+1)*2. +grid[p].r(1,0,2*k-1)*2. + grid[p].r(grid[p].N-1,0,2*k+1)*2.+ grid[p].r(grid[p].N-1,0,2*k-1)*2.
                         +grid[p].r(0,1,2*k+1)*2. +grid[p].r(0,1,2*k-1)*2. + grid[p].r(0,grid[p].M-1,2*k+1)*2.+ grid[p].r(0,grid[p].M-1,2*k-1)*2.
                         +grid[p].r(1,1,2*k+1) +grid[p].r(grid[p].N-1,1,2*k+1)+grid[p].r(1,grid[p].M-1,2*k+1)+grid[p].r(1,1,2*k-1)+grid[p].r(grid[p].N-1,grid[p].M-1,2*k+1)+grid[p].r(grid[p].N-1,1,2*k-1)+grid[p].r(1,grid[p].M-1,2*k-1)+grid[p].r(grid[p].N-1,grid[p].M-1,2*k-1))/64.;

    }
}

void MG3d::CoarseToFine(const int & p)
{
    for(int i=0;i<grid[p].N-1;i++)
    {
        for(int j=0;j<grid[p].M-1;j++)
        {
            for(int k=0;k<grid[p].L-1;k++)
            {
                grid[p+1].u(i*2,j*2,k*2)+=grid[p].u(i,j,k);
                grid[p+1].u(i*2+1,j*2,k*2)+= (grid[p].u(i,j,k)+grid[p].u(i+1,j,k))/2.;
                grid[p+1].u(i*2,j*2+1,k*2)+= (grid[p].u(i,j,k)+grid[p].u(i,j+1,k))/2.;
                grid[p+1].u(i*2,j*2,k*2+1)+= (grid[p].u(i,j,k)+grid[p].u(i,j,k+1))/2.;
                grid[p+1].u(i*2+1,j*2+1,k*2)+= (grid[p].u(i,j,k)+grid[p].u(i+1,j,k)+grid[p].u(i,j+1,k)+grid[p].u(i+1,j+1,k))/4.;
                grid[p+1].u(i*2+1,j*2,k*2+1)+= (grid[p].u(i,j,k)+grid[p].u(i+1,j,k)+grid[p].u(i,j,k+1)+grid[p].u(i+1,j,k+1))/4.;
                grid[p+1].u(i*2,j*2+1,k*2+1)+= (grid[p].u(i,j,k)+grid[p].u(i,j+1,k)+grid[p].u(i,j,k+1)+grid[p].u(i,j+1,k+1))/4.;
                grid[p+1].u(i*2+1,j*2+1,k*2+1)+= (grid[p].u(i,j,k)+grid[p].u(i+1,j,k)+grid[p].u(i,j+1,k)+grid[p].u(i,j,k+1)+
                                            grid[p].u(i+1,j+1,k)+grid[p].u(i+1,j,k+1)+grid[p].u(i,j+1,k+1)+grid[p].u(i+1,j+1,k+1))/8.;
            }
            grid[p+1].u(i*2,j*2,grid[p+1].L-2)+=grid[p].u(i,j,grid[p].L-1);
            grid[p+1].u(i*2+1,j*2,grid[p+1].L-2)+= (grid[p].u(i,j,grid[p].L-1)+grid[p].u(i+1,j,grid[p].L-1))/2.;
            grid[p+1].u(i*2,j*2+1,grid[p+1].L-2)+=(grid[p].u(i,j,grid[p].L-1)+grid[p].u(i,j+1,grid[p].L-1))/2.;
            grid[p+1].u(i*2,j*2,grid[p+1].L-1)+= (grid[p].u(i,j,grid[p].L-1)+grid[p].u(i,j,0))/2.;
            grid[p+1].u(i*2+1,j*2+1,grid[p+1].L-2)+= (grid[p].u(i,j,grid[p].L-1)+grid[p].u(i+1,j,grid[p].L-1)+grid[p].u(i,j+1,grid[p].L-1)+grid[p].u(i+1,j+1,grid[p].L-1))/4.;
            grid[p+1].u(i*2+1,j*2,grid[p+1].L-1)+= (grid[p].u(i,j,grid[p].L-1)+grid[p].u(i+1,j,grid[p].L-1)+grid[p].u(i,j,0)+grid[p].u(i+1,j,0))/4.;
            grid[p+1].u(i*2,j*2+1,grid[p+1].L-1)+= (grid[p].u(i,j,grid[p].L-1)+grid[p].u(i,j+1,grid[p].L-1)+grid[p].u(i,j,0)+grid[p].u(i,j+1,0))/4.;
            grid[p+1].u(i*2+1,j*2+1,grid[p+1].L-1)+= (grid[p].u(i,j,grid[p].L-1)+grid[p].u(i+1,j,grid[p].L-1)+grid[p].u(i,j+1,grid[p].L-1)+grid[p].u(i,j,0)
                                        +grid[p].u(i+1,j+1,grid[p].L-1)+grid[p].u(i+1,j,0)+grid[p].u(i,j+1,0)+grid[p].u(i+1,j+1,0))/8.;
        }
        grid[p+1].u(i*2,grid[p+1].M-2,grid[p+1].L-2)+=grid[p].u(i,grid[p].M-1,grid[p].L-1);
        grid[p+1].u(i*2+1,grid[p+1].M-2,grid[p+1].L-2)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i+1,grid[p].M-1,grid[p].L-1))/2.;
        grid[p+1].u(i*2,grid[p+1].M-1,grid[p+1].L-2)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i,0,grid[p].L-1))/2.;
        grid[p+1].u(i*2,grid[p+1].M-2,grid[p+1].L-1)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i,grid[p].M-1,0))/2.;
        grid[p+1].u(i*2+1,grid[p+1].M-1,grid[p+1].L-2)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i+1,grid[p].M-1,grid[p].L-1)+grid[p].u(i,0,grid[p].L-1)+grid[p].u(i+1,0,grid[p].L-1))/4.;
        grid[p+1].u(i*2+1,grid[p+1].M-2,grid[p+1].L-1)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i+1,grid[p].M-1,grid[p].L-1)+grid[p].u(i,grid[p].M-1,0)+grid[p].u(i+1,grid[p].M-1,0))/4.;
        grid[p+1].u(i*2,grid[p+1].M-1,grid[p+1].L-1)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i,0,grid[p].L-1)+grid[p].u(i,grid[p].M-1,0)+grid[p].u(i,0,0))/4.;
        grid[p+1].u(i*2+1,grid[p+1].M-1,grid[p+1].L-1)+= (grid[p].u(i,grid[p].M-1,grid[p].L-1)+grid[p].u(i+1,grid[p].M-1,grid[p].L-1)+grid[p].u(i,0,grid[p].L-1)+grid[p].u(i,grid[p].M-1,0)
                                    +grid[p].u(i+1,0,grid[p].L-1)+grid[p].u(i+1,grid[p].M-1,0)+grid[p].u(i,0,0)+grid[p].u(i+1,0,0))/8.;
    }
    for(int j=0;j<grid[p].M-1;j++)
    {
        for(int k=0;k<grid[p].L-1;k++)
        {
            grid[p+1].u(grid[p+1].N-2,j*2,k*2)+=grid[p].u(grid[p].N-1,j,k);
            grid[p+1].u(grid[p+1].N-1,j*2,k*2)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(0,j,k))/2.;
            grid[p+1].u(grid[p+1].N-2,j*2+1,k*2)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(grid[p].N-1,j+1,k))/2.;
            grid[p+1].u(grid[p+1].N-2,j*2,k*2+1)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(grid[p].N-1,j,k+1))/2.;
            grid[p+1].u(grid[p+1].N-1,j*2+1,k*2)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(0,j,k)+grid[p].u(grid[p].N-1,j+1,k)+grid[p].u(0,j+1,k))/4.;
            grid[p+1].u(grid[p+1].N-1,j*2,k*2+1)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(0,j,k)+grid[p].u(grid[p].N-1,j,k+1)+grid[p].u(0,j,k+1))/4.;
            grid[p+1].u(grid[p+1].N-2,j*2+1,k*2+1)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(grid[p].N-1,j+1,k)+grid[p].u(grid[p].N-1,j,k+1)+grid[p].u(grid[p].N-1,j+1,k+1))/4.;
            grid[p+1].u(grid[p+1].N-1,j*2+1,k*2+1)+= (grid[p].u(grid[p].N-1,j,k)+grid[p].u(0,j,k)+grid[p].u(grid[p].N-1,j+1,k)+grid[p].u(grid[p].N-1,j,k+1)+
                                        grid[p].u(0,j+1,k)+grid[p].u(0,j,k+1)+grid[p].u(grid[p].N-1,j+1,k+1)+grid[p].u(0,j+1,k+1))/8.;
        }
        grid[p+1].u(grid[p+1].N-2,j*2,grid[p+1].L-2)+=grid[p].u(grid[p].N-1,j,grid[p].L-1);
        grid[p+1].u(grid[p+1].N-1,j*2,grid[p+1].L-2)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(0,j,grid[p].L-1))/2.;
        grid[p+1].u(grid[p+1].N-2,j*2+1,grid[p+1].L-2)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(grid[p].N-1,j+1,grid[p].L-1))/2.;
        grid[p+1].u(grid[p+1].N-2,j*2,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(grid[p].N-1,j,0))/2.;
        grid[p+1].u(grid[p+1].N-1,j*2+1,grid[p+1].L-2)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(0,j,grid[p].L-1)+grid[p].u(grid[p].N-1,j+1,grid[p].L-1)+grid[p].u(0,j+1,grid[p].L-1))/4.;
        grid[p+1].u(grid[p+1].N-1,j*2,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(0,j,grid[p].L-1)+grid[p].u(grid[p].N-1,j,0)+grid[p].u(0,j,0))/4.;
        grid[p+1].u(grid[p+1].N-2,j*2+1,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(grid[p].N-1,j+1,grid[p].L-1)+grid[p].u(grid[p].N-1,j,0)+grid[p].u(grid[p].N-1,j+1,0))/4.;
        grid[p+1].u(grid[p+1].N-1,j*2+1,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,j,grid[p].L-1)+grid[p].u(0,j,grid[p].L-1)+grid[p].u(grid[p].N-1,j+1,grid[p].L-1)+grid[p].u(grid[p].N-1,j,0)+
                                    grid[p].u(0,j+1,grid[p].L-1)+grid[p].u(0,j,0)+grid[p].u(grid[p].N-1,j+1,0)+grid[p].u(0,j+1,0))/8.;

    }
    for(int k=0;k<grid[p].L-1;k++)
    {
            grid[p+1].u(grid[p+1].N-2,grid[p+1].M-2,k*2)+=grid[p].u(grid[p].N-1,grid[p].M-1,k);
            grid[p+1].u(grid[p+1].N-1,grid[p+1].M-2,k*2)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(0,grid[p].M-1,k))/2.;
            grid[p+1].u(grid[p+1].N-2,grid[p+1].M-1,k*2)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(grid[p].N-1,0,k))/2.;
            grid[p+1].u(grid[p+1].N-2,grid[p+1].M-2,k*2+1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(grid[p].N-1,grid[p].M-1,k+1))/2.;
            grid[p+1].u(grid[p+1].N-1,grid[p+1].M-1,k*2)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(0,grid[p].M-1,k)+grid[p].u(grid[p].N-1,0,k)+grid[p].u(0,0,k))/4.;
            grid[p+1].u(grid[p+1].N-1,grid[p+1].M-2,k*2+1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(0,grid[p].M-1,k)+grid[p].u(grid[p].N-1,grid[p].M-1,k+1)+grid[p].u(0,grid[p].M-1,k+1))/4.;
            grid[p+1].u(grid[p+1].N-2,grid[p+1].M-1,k*2+1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(grid[p].N-1,0,k)+grid[p].u(grid[p].N-1,grid[p].M-1,k+1)+grid[p].u(grid[p].N-1,0,k+1))/4.;
            grid[p+1].u(grid[p+1].N-1,grid[p+1].M-1,k*2+1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,k)+grid[p].u(0,grid[p].M-1,k)+grid[p].u(grid[p].N-1,0,k)+grid[p].u(grid[p].N-1,grid[p].M-1,k+1)+
                                        grid[p].u(0,0,k)+grid[p].u(0,grid[p].M-1,k+1)+grid[p].u(grid[p].N-1,0,k+1)+grid[p].u(0,0,k+1))/8.;

    }
    for(int i=0;i<grid[p].N-1;i++)
    {
        for(int k=0;k<grid[p].L-1;k++)
        {
                grid[p+1].u(i*2,grid[p+1].M-2,k*2)+=grid[p].u(i,grid[p].M-1,k);
                grid[p+1].u(i*2+1,grid[p+1].M-2,k*2)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i+1,grid[p].M-1,k))/2.;
                grid[p+1].u(i*2,grid[p+1].M-1,k*2)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i,0,k))/2.;
                grid[p+1].u(i*2,grid[p+1].M-2,k*2+1)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i,grid[p].M-1,k+1))/2.;
                grid[p+1].u(i*2+1,grid[p+1].M-1,k*2)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i+1,grid[p].M-1,k)+grid[p].u(i,0,k)+grid[p].u(i+1,0,k))/4.;
                grid[p+1].u(i*2+1,grid[p+1].M-2,k*2+1)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i+1,grid[p].M-1,k)+grid[p].u(i,grid[p].M-1,k+1)+grid[p].u(i+1,grid[p].M-1,k+1))/4.;
                grid[p+1].u(i*2,grid[p+1].M-1,k*2+1)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i,0,k)+grid[p].u(i,grid[p].M-1,k+1)+grid[p].u(i,0,k+1))/4.;
                grid[p+1].u(i*2+1,grid[p+1].M-1,k*2+1)+= (grid[p].u(i,grid[p].M-1,k)+grid[p].u(i+1,grid[p].M-1,k)+grid[p].u(i,0,k)+grid[p].u(i,grid[p].M-1,k+1)+
                                            grid[p].u(i+1,0,k)+grid[p].u(i+1,grid[p].M-1,k+1)+grid[p].u(i,0,k+1)+grid[p].u(i+1,0,k+1))/8.;

        }
    }
    grid[p+1].u(grid[p+1].N-2,grid[p+1].M-2,grid[p+1].L-2)+=grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1);
    grid[p+1].u(grid[p+1].N-1,grid[p+1].M-2,grid[p+1].L-2)+= (grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(0,grid[p].M-1,grid[p].L-1))/2.;
    grid[p+1].u(grid[p+1].N-2,grid[p+1].M-1,grid[p+1].L-2)+=(grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(grid[p].N-1,0,grid[p].L-1))/2.;
    grid[p+1].u(grid[p+1].N-2,grid[p+1].M-2,grid[p+1].L-1)+=(grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(grid[p].N-1,grid[p].M-1,0))/2.;
    grid[p+1].u(grid[p+1].N-1,grid[p+1].M-1,grid[p+1].L-2)+= (grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(0,grid[p].M-1,grid[p].L-1)+grid[p].u(grid[p].N-1,0,grid[p].L-1)+grid[p].u(0,0,grid[p].L-1))/4.;
    grid[p+1].u(grid[p+1].N-1,grid[p+1].M-2,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(0,grid[p].M-1,grid[p].L-1)+grid[p].u(grid[p].N-1,grid[p].M-1,0)+grid[p].u(0,grid[p].M-1,0))/4.;
    grid[p+1].u(grid[p+1].N-2,grid[p+1].M-1,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(grid[p].N-1,0,grid[p].L-1)+grid[p].u(grid[p].N-1,grid[p].M-1,0)+grid[p].u(grid[p].N-1,0,0))/4.;
    grid[p+1].u(grid[p+1].N-1,grid[p+1].M-1,grid[p+1].L-1)+= (grid[p].u(grid[p].N-1,grid[p].M-1,grid[p].L-1)+grid[p].u(0,grid[p].M-1,grid[p].L-1)+grid[p].u(grid[p].N-1,0,grid[p].L-1)+grid[p].u(grid[p].N-1,grid[p].M-1,0)
                                +grid[p].u(0,0,grid[p].L-1)+grid[p].u(0,grid[p].M-1,0)+grid[p].u(grid[p].N-1,0,0)+grid[p].u(0,0,0))/8.;
}

void MG3d::Solve()
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

void MG3d::Precondition(const Vec3d<double> & b, Vec3d<double> &x)
{
    grid[_size-1].f=b;
    //grid[_size-1].PreSmooth();
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
