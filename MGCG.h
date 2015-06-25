#ifndef MGCG_H
#define MGCG_H

#include "CCF.h"
#include "vector2d.h"

template<typename T>
double GaussSeidel(const CCF<T> & mat, vec<T> & unknown, vec<T> & unknown_new, vec<T> & rhs, const int & ITR_MAX)
{
    double err=1.;
    for(int itr=0;itr<ITR_MAX&&err>1.0E-7;itr++)
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
        cout<<err<<endl;
    }
    return err;
}


//********M: grid number*****************
template<typename T, int M>
class MGCG2d
{
public:
    MGCG2d();
    void Solve();
    virtual void SetMat()=0;
    virtual void SetRHS()=0;
//protected:
    CCF<T> * _mat;
    Vec2d<T> * _rhs;
    Vec2d<T> * _unknown;
    Vec2d<T> * _unknown_new;
};

template<typename T, int M>
MGCG2d<T,M>::MGCG2d()
{
    _mat = new CCF<T> [M];
    _rhs = new Vec2d<T> [M];
    _unknown = new Vec2d<T> [M];
    _unknown_new = new Vec2d<T> [M];
}

template<typename T, int M>
void MGCG2d<T,M>::Solve()
{
    cout<<GaussSeidel(_mat[M-1],_unknown[M-1],_unknown_new[M-1],_rhs[M-1],100)<<endl;
}
#endif // MGCG_H
