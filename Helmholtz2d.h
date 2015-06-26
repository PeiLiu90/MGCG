#ifndef HELMHOLTZ2D_H
#define HELMHOLTZ2D_H

#include "vector2d.h"
#include "MGCG.h"
#include "math.h"
#include "Equation.h"
#include "tools.h"
/*Solve Equation \nabla \eps \nabla \phi + k^2 \phi = \rho
both \eps and k are scalar functions of coordinates*/

//***********************************************************
/*
Nx : mesh number of x direction
Ny : mesh number of y direction
M  : grid number of the coarse and fine grids
O : the coarsest
M-1: the finest
*/
//************************************************************
template<typename T, int Nx, int Ny, int M>
class Helmholtz2d : public MGCG<Vec2d<T>,T,M>
{
public:
    Helmholtz2d(const double & h);
    void SetMat();
    void SetRHS();
    void FineToCoarse(const int & k);
    void CoarseToFine(const int & k);
    T ComputeErr(const int & k) const ;

    int Index(const int & , const int & , const int &);
//private:
    int _size[M][2];
    Vec2d<T> _exact;
    Vec2d<T> * _eps;//this is \eps
    Vec2d<T> * _k;//this is k^2
    double * _x;
    double * _y;
    double * _h;
};

template<typename T, int Nx, int Ny, int M>
Helmholtz2d<T,Nx,Ny,M>::Helmholtz2d(const double & h)
{
    _eps = new Vec2d<T> [M];
    _k = new Vec2d<T> [M];
    _h = new double [M];
    for(int i=0;i<M;i++)
    {
        _h[i] = h * pow(0.5,i);
        _size[i][0]=Nx*pow(2,i);
        _size[i][1]=Ny*pow(2,i);
        this->_mat[i].SetSize(_size[i][0]*_size[i][1]);
        this->_mat[i].SetLength(_size[i][0]*_size[i][1]*3);
        this->_err[i].SetSize(_size[i][0],_size[i][1]);
        this->_rhs[i].SetSize(_size[i][0],_size[i][1]);
        this->_unknown[i].SetSize(_size[i][0],_size[i][1]);
        this->_unknown_new[i].SetSize(_size[i][0],_size[i][1]);
        _eps[i].SetSize(_size[i][0],_size[i][1]);
        _k[i].SetSize(_size[i][0],_size[i][1]);
        for(int j=0;j<_size[i][0];j++)
        {
            for(int k=0;k<_size[i][1];k++)
            {
                _eps[i](j,k)=1.;
                _k[i](j,k)=1.;
            }
        }
    }
    this->_b.SetSize(_size[M-1][0],_size[M-1][1]);
    this->_u.SetSize(_size[M-1][0],_size[M-1][1]);
    this->_r.SetSize(_size[M-1][0],_size[M-1][1]);
    this->_z.SetSize(_size[M-1][0],_size[M-1][1]);
    this->_p.SetSize(_size[M-1][0],_size[M-1][1]);
    this->_Ap.SetSize(_size[M-1][0],_size[M-1][1]);
    _x = new double [_size[M-1][0]];
    _y = new double [_size[M-1][1]];
    for(int i=0;i<_size[M-1][0];i++)
    {
        _x[i]= double(i) * _h[M-1];
    }
    for(int i=0;i<_size[M-1][1];i++)
    {
        _y[i]= double(i) * _h[M-1];
    }
}

template<typename T, int Nx, int Ny, int M>
int Helmholtz2d<T,Nx,Ny,M>::Index(const int & k, const int & i, const int & j)
{
    return mod(i,_size[k][1])*_size[k][0]+mod(j,_size[k][0]);
}

template<typename T, int Nx, int Ny, int M>
void Helmholtz2d<T,Nx,Ny,M>::SetMat()
{
    int index;
    int row;
    int col;
    for(int k=0;k<M;k++)
    {
        index=0;
        col=0;
        for(int i=0;i<_size[k][1];i++)
        {
            for(int j=0;j<_size[k][0];j++)
            {
                this->_mat[k].colptr(col)=index+1;
                col++;
                this->_mat[k].rowind(index)=col;
                this->_mat[k].data(index)= _k[k](i,j)*_h[k]*_h[k]+ 2.*_eps[k](Index(k,i,j))
                    + (_eps[k](Index(k,i-1,j))+_eps[k](Index(k,i+1,j))+_eps[k](Index(k,i,j-1))+_eps[k](Index(k,i,j+1)))*0.5;
                index++;

                row = Index(k,i,j+1) +1;
                if(row > col)
                {
                    this->_mat[k].rowind(index)=row;
                    this->_mat[k].data(index)= -(_eps[k](Index(k,i,j))+_eps[k](Index(k,i,j+1)))*0.5;
                    index++;
                }
                row = Index(k,i,j-1) +1;
                if(row > col)
                {
                    this->_mat[k].rowind(index)=row;
                    this->_mat[k].data(index)=-(_eps[k](Index(k,i,j))+_eps[k](Index(k,i,j-1)))*0.5;
                    index++;
                }
                row = Index(k,i+1,j) +1;
                if(row > col)
                {
                    this->_mat[k].rowind(index)=row;
                    this->_mat[k].data(index)=-(_eps[k](Index(k,i,j))+_eps[k](Index(k,i+1,j)))*0.5;
                    index++;
                }
                row = Index(k,i-1,j) +1;
                if(row > col)
                {
                    this->_mat[k].rowind(index)=row;
                    this->_mat[k].data(index)=-(_eps[k](Index(k,i,j))+_eps[k](Index(k,i-1,j)))*0.5;
                    index++;
                }
            }
        }
        this->_mat[k].colptr(col)=index+1;
    }
};

template<typename T, int Nx, int Ny, int M>
void Helmholtz2d<T,Nx,Ny,M>::SetRHS()
{
    _exact.SetSize(_size[M-1][0],_size[M-1][1]);
    for(int i=0;i<_size[M-1][0];i++)
    {
        for(int j=0;j<_size[M-1][1];j++)
        {
            this->_b(i,j)= source(_x[i],_y[j])*_h[M-1]*_h[M-1];
            _exact(i,j) = exact(_x[i],_y[j]);
        }
    }
}

template<typename T, int Nx, int Ny, int M>
void Helmholtz2d<T,Nx,Ny,M>::FineToCoarse(const int & k)
{
    for(int i=0;i<_size[k-1][0];i++)
    {
        for( int j=0;j<_size[k-1][1];j++)
        {
            this->_rhs[k-1](i,j)=this->_err[k](i*2,j*2);
        }
    }
}

template<typename T, int Nx, int Ny, int M>
void Helmholtz2d<T,Nx,Ny,M>::CoarseToFine(const int & k)
{
    for(int i=0;i<_size[k][0];i++)
    {
        for( int j=0;j<_size[k][1];j++)
        {
            this->_unknown[k+1](i*2,j*2)+=this->_unknown[k](i,j);
            this->_unknown[k+1](i*2+1,j*2)+=(this->_unknown[k](i,j)+this->_unknown[k](i+1,j))*0.5;
            this->_unknown[k+1](i*2,j*2+1)+=(this->_unknown[k](i,j)+this->_unknown[k](i,j+1))*0.5;
            this->_unknown[k+1](i*2+1,j*2+1)+=(this->_unknown[k](i,j)+this->_unknown[k](i+1,j)+this->_unknown[k](i,j+1)+this->_unknown[k](i+1,j+1))*0.25;
        }
    }
}
#endif // HELMHOLTZ2D_H
