#ifndef RBSGS_H
#define RBSGS_H

#include "vector2d.h"

/*
Red Black Symmetric Gauss Seidel smoother
for 2d Elliptic Equation: -\nabla \eps \nabla u + k^2 u = f
both \eps and k are scalar functions of coordinates
*/

class RBSGS
{
public:
    RBSGS(){};
    /*Set up the problem*/
    void SetSize(const int &,const int &);
    void SetCoord(const double &, const double &, const double &, const double &);

    void SetEps(double (*)(const double &, const double &));
    void SetK2(double (*)(const double &, const double &));
    void ComputeD();
    void SetRHS(double (*)(const double &, const double &));

    void PreSmooth();
    void PostSmooth();
    void ComputeAP(const Vec2d<double> &, Vec2d<double> &);
    void Solve(int &, double &);

    void UpdateUr();
    void UpdateUb();
    void ComputeErr();

    Vec2d<double> eps_x;//\eps at x+1/2
    Vec2d<double> eps_y;//\eps at y+1/2
    Vec2d<double> k2;//k^2

    Vec2d<double> D;

    Vec2d<double> u;//solution
    Vec2d<double> f;//rhs
    Vec2d<double> r;//residual
    //THE GRID is DIVIDED into M*(2*N), N*2 is column (x), M is row (y)
    Vec<double> _x;
    Vec<double> _y;
    double hx;
    double hy;
    int N;
    int M;
};

void RBSGS::SetSize(const int & n, const int & m)
{
    N=n;
    M=m;

    eps_x.SetSize(N,M);
    eps_y.SetSize(N,M);

    k2.SetSize(N,M);

    D.SetSize(N,M);
    f.SetSize(N,M);
    u.SetSize(N,M);
    r.SetSize(N,M);
    u=0.;
    r=0.;
    _x.SetSize(N+1);
    _y.SetSize(M+1);
}

void RBSGS::SetCoord(const double & xmin, const double & xmax, const double & ymin, const double & ymax)
{
    hx=(xmax-xmin)/double(N);
    hy=(ymax-ymin)/double(M);
    for(int i=0;i<N+1;i++)
    {
        _x(i)=xmin+hx*double(i);
    }
    for(int j=0;j<M+1;j++)
    {
        _y(j)=ymin+hy*double(j);
    }
}

void RBSGS::SetEps(double (*g)(const double &, const double &))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            eps_x(i,j)=g((_x(i+1)+_x(i))*0.5,_y(j));
            eps_y(i,j)=g(_x(i),(_y(j)+_y(j+1))*0.5);
        }
    }
}

void RBSGS::SetK2(double (*g)(const double &, const double &))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            k2(i,j)=g(_x(i),_y(j));
        }
    }
}

void RBSGS::SetRHS(double (*g)(const double &, const double &))
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            f(i,j)=g(_x(i),_y(j));
        }
    }
}

void RBSGS::ComputeD()
{
    for(int j=0;j<M;j++)
    {
        D(0,j)=k2(0,j)+(eps_x(N-1,j)+eps_x(0,j))/(hx*hx);
        for(int i=1;i<N;i++)
        {
            D(i,j)=k2(i,j)+(eps_x(i-1,j)+eps_x(i,j))/(hx*hx);
        }
    }
    for(int i=0;i<N;i++)
    {
        D(i,0)+=(eps_y(i,M-1)+eps_y(i,0))/(hy*hy);
        for(int j=1;j<M;j++)
        {
            D(i,j)+=(eps_y(i,j-1)+eps_y(i,j))/(hy*hy);
        }
    }
}

void RBSGS::UpdateUr()
{
    for(int i=0;i<N;i+=2)
    {
        u(i,0)=f(i,0)+(u(i,M-1)*eps_y(i,M-1)+u(i,1)*eps_y(i,0))/(hy*hy);
        for(int j=2;j<M;j+=2)
        {
            u(i,j)=f(i,j)+( u(i,j-1)*eps_y(i,j-1) + u(i,j+1)*eps_y(i,j))/(hy*hy);
            u(i+1,j-1)=f(i+1,j-1)+(u(i+1,j-2)*eps_y(i+1,j-2) + u(i+1,j)*eps_y(i+1,j-1))/(hy*hy);
        }
        u(i+1,M-1)=f(i+1,M-1)+(u(i+1,M-2)*eps_y(i+1,M-2)+u(i+1,0)*eps_y(i+1,M-1))/(hy*hy);
    }
    for(int j=0;j<M;j+=2)
    {
        u(0,j)+= (u(N-1,j)*eps_x(N-1,j)+u(1,j)*eps_x(0,j))/(hx*hx);
        for(int i=2;i<N;i+=2)
        {
            u(i,j)+=(u(i-1,j)*eps_x(i-1,j)+u(i+1,j)*eps_x(i,j))/(hx*hx);
            u(i-1,j+1)+=(u(i-2,j+1)*eps_x(i-2,j+1)+u(i,j+1)*eps_x(i-1,j+1))/(hx*hx);
        }
        u(N-1,j+1)+=(u(N-2,j+1)*eps_x(N-2,j+1)+u(0,j+1)*eps_x(N-1,j+1))/(hx*hx);
    }
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            u(i,j)/=D(i,j);
            u(i+1,j+1)/=D(i+1,j+1);
        }
    }
}

void RBSGS::UpdateUb()
{
    for(int i=0;i<N;i+=2)
    {
        u(i+1,0)=f(i+1,0)+(u(i+1,M-1)*eps_y(i+1,M-1)+u(i+1,1)*eps_y(i+1,0))/(hy*hy);
        for(int j=2;j<M;j+=2)
        {
            u(i+1,j)=f(i+1,j)+( u(i+1,j-1)*eps_y(i+1,j-1)+u(i+1,j+1)*eps_y(i+1,j))/(hy*hy);
            u(i,j-1)=f(i,j-1)+( u(i,j-2)*eps_y(i,j-2)+u(i,j)*eps_y(i,j-1))/(hy*hy);
        }
        u(i,M-1)=f(i,M-1)+( u(i,M-2)*eps_y(i,M-2)+u(i,0)*eps_y(i,M-1))/(hy*hy);
    }
    for(int j=0;j<M;j+=2)
    {
        u(0,j+1)+=( u(N-1,j+1)*eps_x(N-1,j+1)+u(1,j+1)*eps_x(0,j+1))/(hx*hx);
        for(int i=2;i<N;i+=2)
        {
            u(i,j+1)+=( u(i-1,j+1)*eps_x(i-1,j+1)+u(i+1,j+1)*eps_x(i,j+1))/(hx*hx);
            u(i-1,j)+=( u(i-2,j)*eps_x(i-2,j)+ u(i,j)*eps_x(i-1,j))/(hx*hx);
        }
        u(N-1,j)+=( u(N-2,j)*eps_x(N-2,j)+u(0,j)*eps_x(N-1,j))/(hx*hx);
    }
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            u(i,j+1)/=D(i,j+1);
            u(i+1,j)/=D(i+1,j);
        }
    }
}

void RBSGS::ComputeErr()
{
    for(int i=0;i<N;i++)
    {
        r(i,0)=f(i,0)-D(i,0)*u(i,0)+(eps_y(i,M-1)*u(i,M-1)+eps_y(i,0)*u(i,1))/(hy*hy);
        for(int j=1;j<M-1;j++)
        {
            r(i,j)=f(i,j)-D(i,j)*u(i,j)+(eps_y(i,j-1)*u(i,j-1)+eps_y(i,j)*u(i,j+1))/(hy*hy);
        }
        r(i,M-1)=f(i,M-1)-D(i,M-1)*u(i,M-1)+(eps_y(i,M-2)*u(i,M-2)+eps_y(i,M-1)*u(i,0))/(hy*hy);
    }
    for(int j=0;j<M;j++)
    {
        r(0,j)+=(eps_x(N-1,j)*u(N-1,j)+eps_x(0,j)*u(1,j))/(hx*hx);
        for(int i=1;i<N-1;i++)
        {
            r(i,j)+= (eps_x(i-1,j)*u(i-1,j)+eps_x(i,j)*u(i+1,j))/(hx*hx);
        }
        r(N-1,j)+= (eps_x(N-2,j)*u(N-2,j)+eps_x(N-1,j)*u(0,j))/(hx*hx);
    }
}

void RBSGS::ComputeAP(const Vec2d<double> & p,Vec2d<double> & Ap)
{
    for(int i=0;i<N;i++)
    {
        Ap(i,0)=D(i,0)*p(i,0)-(eps_y(i,M-1)*p(i,M-1)+eps_y(i,0)*p(i,1))/(hy*hy);
        for(int j=1;j<M-1;j++)
        {
            Ap(i,j)=D(i,j)*p(i,j)-(eps_y(i,j-1)*p(i,j-1)+eps_y(i,j)*p(i,j+1))/(hy*hy);
        }
        Ap(i,M-1)=D(i,M-1)*p(i,M-1)-(eps_y(i,M-2)*p(i,M-2)+eps_y(i,M-1)*p(i,0))/(hy*hy);
    }
    for(int j=0;j<M;j++)
    {
        Ap(0,j)-=(eps_x(N-1,j)*p(N-1,j)+eps_x(0,j)*p(1,j))/(hx*hx);
        for(int i=1;i<N-1;i++)
        {
            Ap(i,j)-= (eps_x(i-1,j)*p(i-1,j)+eps_x(i,j)*p(i+1,j))/(hx*hx);
        }
        Ap(N-1,j)-= (eps_x(N-2,j)*p(N-2,j)+eps_x(N-1,j)*p(0,j))/(hx*hx);
    }
}

void RBSGS::Solve(int & itr, double & err)
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
        cout<<itr<<"  "<<err<<endl;
    }
}

void RBSGS::PreSmooth()
{
    for(int i=0;i<N;i+=2)
    {
        for(int j=0;j<M;j+=2)
        {
            u(i,j)=f(i,j)/D(i,j);
            u(i+1,j+1)=f(i+1,j+1)/D(i+1,j+1);
        }
    }
    UpdateUb();
    UpdateUr();
    ComputeErr();
}

void RBSGS::PostSmooth()
{
    UpdateUr();
    UpdateUb();
    UpdateUr();
}













#endif // RBSGS_H
