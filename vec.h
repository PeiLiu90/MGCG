#ifndef VECTOR_H
#define VECTOR_H

#include "math.h"
#include <iostream>
using namespace std;

template<typename T>
class vec
{
public:
    virtual T & operator()(const int & ) = 0;
    virtual const T & operator()(const int &) const = 0;
};

template<typename E, typename T>
class VecExpression
{
public:
    T operator()(const int & i) const{return static_cast<const E&>(*this)(i);}
    int size() const {return static_cast<const E &>(*this).size();}

    operator E&() { return static_cast<E&>(*this);}
    operator E const&() const { return static_cast<const E&>(*this);}
};

template<typename T>
class Vec : public VecExpression<Vec<T> , T >, public vec<T>
{
public:
    Vec(){};
    Vec(const Vec<T> & u)
    {
        _size= u.size();
        _data=new T [_size];
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
    }
    template<typename E>
    Vec(const VecExpression<E,T> & u)
    {
        _size= u.size();
        _data=new T [_size];
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
    }
    Vec(const int & n){_size = n; _data = new T [n];}

    void SetSize(const int & n) {_size = n; _data = new T [n];}

    T & operator()(const int &i) {return _data[i];}
    const T & operator()(const int &i) const {return _data[i];}

    const int & size() const {return _size;}

    const Vec<T> & operator=(const Vec<T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
        return *this;
    }
    template<typename E>
    const Vec<T> & operator=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec<T> & operator+=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]+=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec<T> & operator-=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]-=u(i);
        }
        return *this;
    }

    const Vec<T> & operator*=(const T & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]*=u;
        }
        return *this;
    }

    const Vec<T> & operator/=(const T & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]/=u;
        }
        return *this;
    }
private:
    T * _data;
    int _size;
};

template<typename E1, typename E2, typename T>
class VecPlus : public VecExpression<VecPlus<E1,E2,T> ,T>
{
public:
    VecPlus(const E1 & u, const E2 & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)+_v(i);}
    int size() const {return _u.size();}
private:
    const E1 & _u;
    const E2 & _v;
};

template<typename E1, typename E2, typename T>
VecPlus<E1,E2,T> operator+(const VecExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    return VecPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
class VecMinus : public VecExpression<VecMinus<E1,E2,T> ,T>
{
public:
    VecMinus(const E1 & u, const E2 & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)-_v(i);}
    int size() const {return _u.size();}
private:
    const E1 & _u;
    const E2 & _v;
};

template<typename E1, typename E2, typename T>
VecMinus<E1,E2,T> operator-(const VecExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    return VecMinus<E1,E2,T>(u,v);
}

template<typename E, typename T>
class VecPlusScalar : public VecExpression<VecPlusScalar<E,T> ,T>
{
public:
    VecPlusScalar(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)+_v;}
    int size() const {return _u.size();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
VecPlusScalar<E,T> operator+(const VecExpression<E,T> & u, const T & v)
{
    return VecPlusScalar<E,T>(u,v);
}

template<typename E, typename T>
VecPlusScalar<E,T> operator+(const T & v, const VecExpression<E,T> & u)
{
    return VecPlusScalar<E,T>(u,v);
}

template<typename E, typename T>
class VecMinusScalar : public VecExpression<VecMinusScalar<E,T> ,T>
{
public:
    VecMinusScalar(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)-_v;}
    int size() const {return _u.size();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
VecMinusScalar<E,T> operator-(const VecExpression<E,T> & u, const T & v)
{
    return VecMinusScalar<E,T>(u,v);
}

template<typename E, typename T>
class ScalarMinusVec : public VecExpression<ScalarMinusVec<E,T> ,T>
{
public:
    ScalarMinusVec(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _v-_u(i);}
    int size() const {return _u.size();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
ScalarMinusVec<E,T> operator-(const T & v,const VecExpression<E,T> & u)
{
    return ScalarMinusVec<E,T>(u,v);
}

template<typename E1, typename E2, typename T>
T operator*(const VecExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    T tmp=0;
    for(int i=0;i<u.size();i++)
    {
        tmp+=u(i)*v(i);
    }
    return tmp;
}

template<typename E, typename T>
T norm(const VecExpression<E,T> & u)
{
    T tmp=0;
    for(int i=0;i<u.size();i++)
    {
        tmp+=u(i)*u(i);
    }
    tmp/=T(u.size());
    return sqrt(tmp);
}

template<typename E, typename T>
T normmax(const VecExpression<E,T> & u)
{
    T tmp=0;
    for(int i=0;i<u.size();i++)
    {
        tmp= tmp>fabs(u(i))? tmp : fabs(u(i));
    }
    return tmp;
}

template<typename E1, typename E2, typename T>
void product(const VecExpression<E1,T> & u, const VecExpression<E2,T> & v, double & p)
{
    p=0;
    for(int i=0;i<u.size();i++)
    {
        p+=u(i)*v(i);
    }
}

template<typename E, typename T>
ostream & operator<<(ostream & cout, const VecExpression<E,T> & u)
{
    for(int i=0;i<u.size();i++)
    {
        cout<<u(i)<<endl;
    }
    return cout;
}


#endif // VECTOR_H
