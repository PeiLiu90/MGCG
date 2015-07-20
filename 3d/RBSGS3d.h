#ifndef RBSGS3d_H
#define RBSGS3d_H

#include "vector3d.h"

/*
Red Black Symmetric Gauss Seidel smoother
for 3d Elliptic Equation: -\nabla \eps \nabla u + k^2 u = f
both \eps and k are scalar functions of coordinates
*/

class RBSGS3d
{
public:
    RBSGS3d(){};
    /*Set up the problem*/
    void SetSize(const int &,const int &,const int &);
    void SetCoord(const double &, const double &, const double &, const double &, const double &, const double &);

    void SetEps(double (*)(const double &, const double &, const double &));
    void SetK2(double (*)(const double &, const double &, const double &));
    void ComputeD();
    void SetRHS(double (*)(const double &, const double &, const double &));
    void SetEx(double (*)(const double &, const double &, const double &));

    void PreSmooth();
    void PostSmooth();
    void ComputeAP(const Vec3d<double> &, Vec3d<double> &);
    void Solve(int &, double &);

    void UpdateUr();
    void UpdateUb();
    void ComputeErr();

    Vec3d<double> eps_x;//\eps at x+1/2
    Vec3d<double> eps_y;//\eps at y+1/2
    Vec3d<double> eps_z;//\eps at z+1/2
    Vec3d<double> k2;//k^2

    Vec3d<double> D;
    Vec3d<double> ex;
    Vec3d<double> u;//solution
    Vec3d<double> f;//rhs
    Vec3d<double> r;//residual
    //THE GRID is DIVIDED into M*(2*N), N*2 is column (x), M is row (y)
    Vec<double> _x;
    Vec<double> _y;
    Vec<double> _z;
    double hx;
    double hy;
    double hz;
    int N;
    int M;
    int L;
};

void RBSGS3d::SetSize(const int & n, const int & m, const int & l)
{
    N=n;
    M=m;
    L=l;

    eps_x.SetSize(N,M,L);
    eps_y.SetSize(N,M,L);
    eps_z.SetSize(N,M,L);

    k2.SetSize(N,M,L);
    ex.SetSize(N,M,L);
    D.SetSize(N,M,L);
    f.SetSize(N,M,L);
    u.SetSize(N,M,L);
    r.SetSize(N,M,L);
    u=0.;
    r=0.;
    _x.SetSize(N+1);
    _y.SetSize(M+1);
    _z.SetSize(L+1);
}

void RBSGS3d::SetCoord(const double & xmin, const double & xmax, const double & ymin, const double & ymax, const double & zmin, const double & zmax)
{
    hx=(xmax-xmin)/double(N);
    hy=(ymax-ymin)/double(M);
    hz=(zmax-zmin)/double(L);
    for(int i=0;i<=N;i++)
    {
        _x(i)=xmin+hx*double(i);
    }
    for(int j=0;j<=M;j++)
    {
        _y(j)=ymin+hy*double(j);
    }
    for(int k=0;k<=L;k++)
    {
        _z(k)=zmin+hz*double(k);
    }
}

void RBSGS3d::SetEps(double (*g)(const double &, const double &, const double & ))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            for(int k=0;k<L;k++)
            {
                eps_x(i,j,k)=g((_x(i+1)+_x(i))*0.5,_y(j),_z(k));
                eps_y(i,j,k)=g(_x(i),(_y(j)+_y(j+1))*0.5,_z(k));
                eps_z(i,j,k)=g(_x(i),_y(j),(_z(k)+_z(k+1))*0.5);
            }
        }
    }
}

void RBSGS3d::SetK2(double (*g)(const double &, const double &, const double &))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            for(int k=0;k<L;k++)
            {
                k2(i,j,k)=g(_x(i),_y(j),_z(k));
            }
        }
    }
}

void RBSGS3d::SetRHS(double (*g)(const double &, const double &, const double &))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            for(int k=0;k<L;k++)
            {
                f(i,j,k)=g(_x(i),_y(j),_z(k));
            }
        }
    }
}

void RBSGS3d::SetEx(double (*g)(const double &, const double &, const double &))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            for(int k=0;k<L;k++)
            {
                ex(i,j,k)=g(_x(i),_y(j),_z(k));
            }
        }
    }
}

void RBSGS3d::ComputeD()
{
    for(int k=0;k<L;k++)
    {
        for(int j=0;j<M;j++)
        {
            D(0,j,k)=k2(0,j,k)+(eps_x(N-1,j,k)+eps_x(0,j,k))/(hx*hx);
            for(int i=1;i<N;i++)
            {
                D(i,j,k)=k2(i,j,k)+(eps_x(i-1,j,k)+eps_x(i,j,k))/(hx*hx);
            }
        }
    }
    for(int k=0;k<L;k++)
    {
        for(int i=0;i<N;i++)
        {
            D(i,0,k)+=(eps_y(i,M-1,k)+eps_y(i,0,k))/(hy*hy);
            for(int j=1;j<M;j++)
            {
                D(i,j,k)+=(eps_y(i,j-1,k)+eps_y(i,j,k))/(hy*hy);
            }
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            D(i,j,0)+=(eps_z(i,j,L-1)+eps_z(i,j,0))/(hz*hz);
            for(int k=1;k<L;k++)
            {
                D(i,j,k)+=(eps_z(i,j,k-1)+eps_z(i,j,k))/(hz*hz);
            }
        }
    }
}

void RBSGS3d::UpdateUr()
{
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            u(i,j,0)=f(i,j,0) + (u(i,j,L-1)*eps_z(i,j,L-1)+u(i,j,1)*eps_z(i,j,0))/(hz*hz);
            u(i+1,j+1,0)=f(i+1,j+1,0) + (u(i+1,j+1,L-1)*eps_z(i+1,j+1,L-1)+u(i+1,j+1,1)*eps_z(i+1,j+1,0))/(hz*hz);
            for(int k=2;k<L;k+=2)
            {
                u(i,j,k)= f(i,j,k) + (u(i,j,k-1)*eps_z(i,j,k-1)+u(i,j,k+1)*eps_z(i,j,k))/(hz*hz);
                u(i+1,j+1,k)= f(i+1,j+1,k) + (u(i+1,j+1,k-1)*eps_z(i+1,j+1,k-1)+u(i+1,j+1,k+1)*eps_z(i+1,j+1,k))/(hz*hz);
                u(i,j+1,k-1)= f(i,j+1,k-1) + (u(i,j+1,k-2)*eps_z(i,j+1,k-2)+u(i,j+1,k)*eps_z(i,j+1,k-1))/(hz*hz);
                u(i+1,j,k-1)= f(i+1,j,k-1) + (u(i+1,j,k-2)*eps_z(i+1,j,k-2)+u(i+1,j,k)*eps_z(i+1,j,k-1))/(hz*hz);
            }
            u(i,j+1,L-1)=f(i,j+1,L-1) + (u(i,j+1,L-2)*eps_z(i,j+1,L-2)+u(i,j+1,0)*eps_z(i,j+1,L-1))/(hz*hz);
            u(i+1,j,L-1)=f(i+1,j,L-1) + (u(i+1,j,L-2)*eps_z(i+1,j,L-2)+u(i+1,j,0)*eps_z(i+1,j,L-1))/(hz*hz);
        }
    }
    for(int i=0;i<N;i+=2)
    {
        for(int k=0;k<L;k+=2)
        {
            u(i,0,k) += (u(i,M-1,k)*eps_y(i,M-1,k)+u(i,1,k)*eps_y(i,0,k))/(hy*hy);
            u(i+1,0,k+1) += (u(i+1,M-1,k+1)*eps_y(i+1,M-1,k+1)+u(i+1,1,k+1)*eps_y(i+1,0,k+1))/(hy*hy);
            for(int j=2;j<M;j+=2)
            {
                u(i,j,k) += (u(i,j-1,k)*eps_y(i,j-1,k)+u(i,j+1,k)*eps_y(i,j,k))/(hy*hy);
                u(i+1,j,k+1) += (u(i+1,j-1,k+1)*eps_y(i+1,j-1,k+1)+ u(i+1,j+1,k+1)*eps_y(i+1,j,k+1))/(hy*hy);
                u(i+1,j-1,k) += (u(i+1,j-2,k)*eps_y(i+1,j-2,k)+ u(i+1,j,k)*eps_y(i+1,j-1,k))/(hy*hy);
                u(i,j-1,k+1) += (u(i,j-2,k+1)*eps_y(i,j-2,k+1)+ u(i,j,k+1)*eps_y(i,j-1,k+1))/(hy*hy);
            }
            u(i+1,M-1,k) += (u(i+1,M-2,k)*eps_y(i+1,M-2,k)+u(i+1,0,k)*eps_y(i+1,M-1,k))/(hy*hy);
            u(i,M-1,k+1) += (u(i,M-2,k+1)*eps_y(i,M-2,k+1)+u(i,0,k+1)*eps_y(i,M-1,k+1))/(hy*hy);
        }
    }
    for(int j=0;j<M;j+=2)
    {
        for(int k=0;k<L;k+=2)
        {
            u(0,j,k) += (u(N-1,j,k)*eps_x(N-1,j,k)+u(1,j,k)*eps_x(0,j,k))/(hx*hx);
            u(0,j+1,k+1) += (u(N-1,j+1,k+1)*eps_x(N-1,j+1,k+1)+u(1,j+1,k+1)*eps_x(0,j+1,k+1))/(hx*hx);
            for(int i=2;i<N;i+=2)
            {
                u(i,j,k) += (u(i-1,j,k)*eps_x(i-1,j,k)+u(i+1,j,k)*eps_x(i,j,k))/(hx*hx);
                u(i-1,j+1,k) += (u(i-2,j+1,k)*eps_x(i-2,j+1,k)+u(i,j+1,k)*eps_x(i-1,j+1,k))/(hx*hx);
                u(i-1,j,k+1) += (u(i-2,j,k+1)*eps_x(i-2,j,k+1)+u(i,j,k+1)*eps_x(i-1,j,k+1))/(hx*hx);
                u(i,j+1,k+1) += (u(i-1,j+1,k+1)*eps_x(i-1,j+1,k+1)+u(i+1,j+1,k+1)*eps_x(i,j+1,k+1))/(hx*hx);
            }
            u(N-1,j+1,k)+= (u(N-2,j+1,k)*eps_x(N-2,j+1,k)+u(0,j+1,k)*eps_x(N-1,j+1,k))/(hx*hx);
            u(N-1,j,k+1)+= (u(N-2,j,k+1)*eps_x(N-2,j,k+1)+u(0,j,k+1)*eps_x(N-1,j,k+1))/(hx*hx);
        }
    }
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            for(int k=0;k<L;k+=2)
            {
                u(i,j,k)/=D(i,j,k);
                u(i+1,j+1,k)/=D(i+1,j+1,k);
                u(i+1,j,k+1)/=D(i+1,j,k+1);
                u(i,j+1,k+1)/=D(i,j+1,k+1);
            }
        }
    }
}

void RBSGS3d::UpdateUb()
{
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            u(i,j+1,0)=f(i,j+1,0) + (u(i,j+1,L-1)*eps_z(i,j+1,L-1)+u(i,j+1,1)*eps_z(i,j+1,0))/(hz*hz);
            u(i+1,j,0)=f(i+1,j,0) + (u(i+1,j,L-1)*eps_z(i+1,j,L-1)+u(i+1,j,1)*eps_z(i+1,j,0))/(hz*hz);
            for(int k=2;k<L;k+=2)
            {
                u(i,j+1,k)= f(i,j+1,k) + (u(i,j+1,k-1)*eps_z(i,j+1,k-1)+u(i,j+1,k+1)*eps_z(i,j+1,k))/(hz*hz);
                u(i+1,j,k)= f(i+1,j,k) + (u(i+1,j,k-1)*eps_z(i+1,j,k-1)+u(i+1,j,k+1)*eps_z(i+1,j,k))/(hz*hz);
                u(i,j,k-1)= f(i,j,k-1) + (u(i,j,k-2)*eps_z(i,j,k-2)+u(i,j,k)*eps_z(i,j,k-1))/(hz*hz);
                u(i+1,j+1,k-1)= f(i+1,j+1,k-1) + (u(i+1,j+1,k-2)*eps_z(i+1,j+1,k-2)+u(i+1,j+1,k)*eps_z(i+1,j+1,k-1))/(hz*hz);
            }
            u(i+1,j+1,L-1)=f(i+1,j+1,L-1) + (u(i+1,j+1,L-2)*eps_z(i+1,j+1,L-2)+u(i+1,j+1,0)*eps_z(i+1,j+1,L-1))/(hz*hz);
            u(i,j,L-1)=f(i,j,L-1) + (u(i,j,L-2)*eps_z(i,j,L-2)+u(i,j,0)*eps_z(i,j,L-1))/(hz*hz);
        }
    }
    for(int i=0;i<N;i+=2)
    {
        for(int k=0;k<L;k+=2)
        {
            u(i,0,k+1) += (u(i,M-1,k+1)*eps_y(i,M-1,k+1)+u(i,1,k+1)*eps_y(i,0,k+1))/(hy*hy);
            u(i+1,0,k) += (u(i+1,M-1,k)*eps_y(i+1,M-1,k)+u(i+1,1,k)*eps_y(i+1,0,k))/(hy*hy);
            for(int j=2;j<M;j+=2)
            {
                u(i,j-1,k) += (u(i,j-2,k)*eps_y(i,j-2,k)+u(i,j,k)*eps_y(i,j-1,k))/(hy*hy);
                u(i+1,j-1,k+1) += (u(i+1,j-2,k+1)*eps_y(i+1,j-2,k+1)+ u(i+1,j,k+1)*eps_y(i+1,j-1,k+1))/(hy*hy);
                u(i+1,j,k) += (u(i+1,j-1,k)*eps_y(i+1,j-1,k)+ u(i+1,j+1,k)*eps_y(i+1,j,k))/(hy*hy);
                u(i,j,k+1) += (u(i,j-1,k+1)*eps_y(i,j-1,k+1)+ u(i,j+1,k+1)*eps_y(i,j,k+1))/(hy*hy);
            }
            u(i,M-1,k) += (u(i,M-2,k)*eps_y(i,M-2,k)+u(i,0,k)*eps_y(i,M-1,k))/(hy*hy);
            u(i+1,M-1,k+1) += (u(i+1,M-2,k+1)*eps_y(i+1,M-2,k+1)+u(i+1,0,k+1)*eps_y(i+1,M-1,k+1))/(hy*hy);
        }
    }
    for(int j=0;j<M;j+=2)
    {
        for(int k=0;k<L;k+=2)
        {
            u(0,j,k+1) += (u(N-1,j,k+1)*eps_x(N-1,j,k+1)+u(1,j,k+1)*eps_x(0,j,k+1))/(hx*hx);
            u(0,j+1,k) += (u(N-1,j+1,k)*eps_x(N-1,j+1,k)+u(1,j+1,k)*eps_x(0,j+1,k))/(hx*hx);
            for(int i=2;i<N;i+=2)
            {
                u(i-1,j,k) += (u(i-2,j,k)*eps_x(i-2,j,k)+u(i,j,k)*eps_x(i-1,j,k))/(hx*hx);
                u(i,j+1,k) += (u(i-1,j+1,k)*eps_x(i-1,j+1,k)+u(i+1,j+1,k)*eps_x(i,j+1,k))/(hx*hx);
                u(i,j,k+1) += (u(i-1,j,k+1)*eps_x(i-1,j,k+1)+u(i+1,j,k+1)*eps_x(i,j,k+1))/(hx*hx);
                u(i-1,j+1,k+1) += (u(i-2,j+1,k+1)*eps_x(i-2,j+1,k+1)+u(i,j+1,k+1)*eps_x(i-1,j+1,k+1))/(hx*hx);
            }
            u(N-1,j,k)+= (u(N-2,j,k)*eps_x(N-2,j,k)+u(0,j,k)*eps_x(N-1,j,k))/(hx*hx);
            u(N-1,j+1,k+1)+= (u(N-2,j+1,k+1)*eps_x(N-2,j+1,k+1)+u(0,j+1,k+1)*eps_x(N-1,j+1,k+1))/(hx*hx);
        }
    }
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            for(int k=0;k<L;k+=2)
            {
                u(i,j,k+1)/=D(i,j,k+1);
                u(i+1,j+1,k+1)/=D(i+1,j+1,k+1);
                u(i+1,j,k)/=D(i+1,j,k);
                u(i,j+1,k)/=D(i,j+1,k);
            }
        }
    }
}


void RBSGS3d::ComputeErr()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            r(i,j,0)=f(i,j,0)-D(i,j,0)*u(i,j,0)+(eps_z(i,j,M-1)*u(i,j,M-1)+eps_z(i,j,0)*u(i,j,1))/(hz*hz);
            for(int k=1;k<L-1;k++)
            {
                r(i,j,k)=f(i,j,k)-D(i,j,k)*u(i,j,k)+(eps_z(i,j,k-1)*u(i,j,k-1)+eps_z(i,j,k)*u(i,j,k+1))/(hz*hz);
            }
            r(i,j,L-1)=f(i,j,L-1)-D(i,j,L-1)*u(i,j,L-1)+(eps_z(i,j,L-2)*u(i,j,L-2)+eps_z(i,j,L-1)*u(i,j,0))/(hz*hz);
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int k=0;k<L;k++)
        {
            r(i,0,k)+= (eps_y(i,M-1,k)*u(i,M-1,k)+eps_y(i,0,k)*u(i,1,k))/(hy*hy);
            for(int j=1;j<M-1;j++)
            {
                r(i,j,k)+= (eps_y(i,j-1,k)*u(i,j-1,k)+eps_y(i,j,k)*u(i,j+1,k))/(hy*hy);
            }
            r(i,M-1,k)+= (eps_y(i,M-2,k)*u(i,M-2,k)+eps_y(i,M-1,k)*u(i,0,k))/(hy*hy);
        }
    }
    for(int j=0;j<M;j++)
    {
        for(int k=0;k<L;k++)
        {
            r(0,j,k)+= (eps_x(N-1,j,k)*u(N-1,j,k)+eps_x(0,j,k)*u(1,j,k))/(hx*hx);
            for(int i=1;i<N-1;i++)
            {
                r(i,j,k)+= (eps_x(i-1,j,k)*u(i-1,j,k)+eps_x(i,j,k)*u(i+1,j,k))/(hx*hx);
            }
            r(N-1,j,k)+= (eps_x(N-2,j,k)*u(N-2,j,k)+eps_x(N-1,j,k)*u(0,j,k))/(hx*hx);
        }
    }
}

void RBSGS3d::ComputeAP(const Vec3d<double> & p,Vec3d<double> & Ap)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            Ap(i,j,0)=D(i,j,0)*p(i,j,0)-(eps_z(i,j,M-1)*p(i,j,M-1)+eps_z(i,j,0)*p(i,j,1))/(hz*hz);
            for(int k=1;k<L-1;k++)
            {
                Ap(i,j,k)=D(i,j,k)*p(i,j,k)-(eps_z(i,j,k-1)*p(i,j,k-1)+eps_z(i,j,k)*p(i,j,k+1))/(hz*hz);
            }
            Ap(i,j,L-1)=D(i,j,L-1)*p(i,j,L-1)-(eps_z(i,j,L-2)*p(i,j,L-2)+eps_z(i,j,L-1)*p(i,j,0))/(hz*hz);
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int k=0;k<L;k++)
        {
            Ap(i,0,k)-= (eps_y(i,M-1,k)*p(i,M-1,k)+eps_y(i,0,k)*p(i,1,k))/(hy*hy);
            for(int j=1;j<M-1;j++)
            {
                Ap(i,j,k)-= (eps_y(i,j-1,k)*p(i,j-1,k)+eps_y(i,j,k)*p(i,j+1,k))/(hy*hy);
            }
            Ap(i,M-1,k)-= (eps_y(i,M-2,k)*p(i,M-2,k)+eps_y(i,M-1,k)*p(i,0,k))/(hy*hy);
        }
    }
    for(int j=0;j<M;j++)
    {
        for(int k=0;k<L;k++)
        {
            Ap(0,j,k)-= (eps_x(N-1,j,k)*p(N-1,j,k)+eps_x(0,j,k)*p(1,j,k))/(hx*hx);
            for(int i=1;i<N-1;i++)
            {
                Ap(i,j,k)-= (eps_x(i-1,j,k)*p(i-1,j,k)+eps_x(i,j,k)*p(i+1,j,k))/(hx*hx);
            }
            Ap(N-1,j,k)-= (eps_x(N-2,j,k)*p(N-2,j,k)+eps_x(N-1,j,k)*p(0,j,k))/(hx*hx);
        }
    }
}

void RBSGS3d::Solve(int & itr, double & err)
{
    itr=0;
    err=1.;
    while(err>1.E-5)
    {
        UpdateUr();
        UpdateUb();
        ComputeErr();
        itr++;
        err=norm(r);
        //cout<<itr<<"  "<<err<<endl;
    }
}

void RBSGS3d::PreSmooth()
{
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            for(int k=0;k<L;k+=2)
            {
                u(i,j,k)=f(i,j,k)/D(i,j,k);
                u(i,j+1,k+1)=f(i,j+1,k+1)/D(i,j+1,k+1);
                u(i+1,j+1,k)=f(i+1,j+1,k)/D(i+1,j+1,k);
                u(i+1,j,k+1)=f(i+1,j,k+1)/D(i+1,j,k+1);
            }
        }
    }
    UpdateUb();
    UpdateUr();
    ComputeErr();
}

void RBSGS3d::PostSmooth()
{
    UpdateUr();
    UpdateUb();
    UpdateUr();
}

#endif // RBSGS3d_H
