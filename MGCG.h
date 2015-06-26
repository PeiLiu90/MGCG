#ifndef MGCG_H
#define MGCG_H

#include "CCF.h"

const int GS_ITR_MAX=100;

template<typename E, typename T>
double GaussSeidel(const CCF<T> & mat, E & unknown, E & unknown_new, const E & rhs, const int & ITR_MAX)
{
    double err=1.;
    for(int itr=0;itr<ITR_MAX&&err>1.0E-8;itr++)
    {
        err=0.;
        for(int i=0;i<mat.size();i++)
        {
            unknown_new(i)=rhs(i);
        }
        for(int i=0;i<mat.size();i++)
        {
            //Diagonal is located at colptr[i]-1
            for(int j=mat.colptr(i);j<mat.colptr(i+1)-1;j++)
            {
                unknown_new(i) -= mat.data(j)* unknown( mat.rowind(j)-1 );
            }
            unknown_new(i)/=mat.data(mat.colptr(i)-1);
            for(int j=mat.colptr(i);j<mat.colptr(i+1)-1;j++)
            {
                unknown_new( mat.rowind(j)-1 ) -= unknown_new(i) * mat.data(j);
            }
        }
        for(int i=0;i<mat.size();i++)
        {
            err+=(unknown(i)-unknown_new(i))*(unknown(i)-unknown_new(i));
            unknown(i)=unknown_new(i);
        }
        err/=double(mat.size());
        err = sqrt(err);
        //cout<<itr<<"  "<<err<<endl;
    }
    return err;
}

//***********************************************************
/*
Nx : mesh number of x direction
Ny : mesh number of y direction
M  : grid number of the coarse and fine grids
O : the coarsest
M-1: the finest
*/
//***************************M: grid number*****************
template<typename E, typename T, int M>
class MGCG
{
public:
    MGCG();
    void Solve();


    void Solve(const int & k);//solve in k-th grid
    void Solve(const E & );//solve for given rhs
    virtual void FineToCoarse(const int & k)=0;//from k to k-1
    virtual void CoarseToFine(const int & k)=0;//from k to k+1
    virtual void SetMat()=0;
    virtual void SetRHS()=0;
//protected:
//Used for Conjugate Gradient
    E _b;
    E _u;
    E _r;
    E _z;
    E _p;
    E _Ap;
//Used for Multi-grid
    CCF<T> * _mat;
    E * _rhs;
    E * _err;
    E * _unknown;
    E * _unknown_new;
};

template<typename E, typename T, int M>
MGCG<E,T,M>::MGCG()
{
    _mat = new CCF<T> [M];
    _err = new E [M];
    _rhs = new E [M];
    _unknown = new E [M];
    _unknown_new = new E [M];
}

template<typename E, typename T, int M>
void MGCG<E,T,M>::Solve(const int & k)
{
    cout<<k<<"  "<<GaussSeidel(_mat[k],_unknown[k],_unknown_new[k],_rhs[k],GS_ITR_MAX)<<endl;
    for(int i=0;i<_mat[k].size();i++)
    {
        _err[k](i) = _rhs[k](i)-_mat[k].data( _mat[k].colptr(i)-1)*_unknown[k](i);
    }
    for(int i=0;i<_mat[k].size();i++)
    {
        //Diagonal is located at colptr[i]-1
        for(int j=_mat[k].colptr(i);j<_mat[k].colptr(i+1)-1;j++)
        {
            _err[k](i) -= _mat[k].data(j)* _unknown[k]( _mat[k].rowind(j)-1 );
            _err[k]( _mat[k].rowind(j)-1 ) -= _unknown[k](i) * _mat[k].data(j);
        }
    }
}

template<typename E, typename T, int M>
void MGCG<E,T,M>::Solve()
{
    Solve(_b);
}

template<typename E, typename T, int M>
void MGCG<E,T,M>::Solve(const E & rhs)
{
    _rhs[M-1]=rhs;
    _unknown[M-1]=_rhs[M-1];
    Solve(M-1);
    for(int itr=0;itr<40;itr++)
    {for(int k=M-1;k>0;k--)
    {
        FineToCoarse(k);
        _unknown[k-1]=_rhs[k-1];
        Solve(k-1);
    }
    for(int k=0;k<M-1;k++)
    {
        CoarseToFine(k);
        Solve(k+1);
    }
    }
}
#endif // MGCG_H
