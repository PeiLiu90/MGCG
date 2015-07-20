#ifndef VECTOR3d_H
#define VECTOR3d_H

#include "math.h"
#include "vec.h"

template<typename E, typename T>
class Vec3dExpression : public VecExpression<E,T>
{
public:
    T operator()(const int & i) const{return static_cast<const E&>(*this)(i);}
    T operator()(const int & i , const int & j, const int & k) const {return static_cast<const E&>(*this)(i,j,k);}
    int size() const {return static_cast<const E &>(*this).size();}
    int size1() const {return static_cast<const E &>(*this).size1();}
    int size2() const {return static_cast<const E &>(*this).size2();}
    int size3() const {return static_cast<const E &>(*this).size3();}

    operator E&() { return static_cast<E&>(*this);}
    operator E const&() const { return static_cast<const E&>(*this);}
};

template<typename T>
class Vec3d : public Vec3dExpression<Vec3d<T>,T>, public vec<T>
{
public:
    Vec3d(){};
    Vec3d(const Vec3d<T> & u)
    {
        _N=u.size1();
        _M=u.size2();
        _L=u.size3();
        _size= _N*_M*_L;
        _data = new T [_size];
        for(int i=0;i<_N*_M*_L;i++)
        {
            _data[i]=u(i);
        }
    }
    template<typename E>
    Vec3d(const Vec3dExpression<E,T> & u)
    {
        _N= u.size1();
        _M= u.size2();
        _L= u.size3();
        _size= _N*_M*_L;
        _data = new T [_size];
        for(int i=0;i<_N*_M*_L;i++)
        {
            _data[i]=u(i);
        }
    }
    Vec3d(const int & n, const int & m, const int & l) {_N=n;_M=m;_L=l;_size= _N*_M*_L; _data = new T [_size];}
    void SetSize(const int & N, const int & M, const int & L){_N= N;_M= M;_L=L;_size= _N*_M*_L; _data = new T [_size];}

    T & operator()(const int &i) {return _data[i];}
    const T & operator()(const int &i) const {return _data[i];}

    T & operator() (const int & i, const int & j, const int & k) {return _data[(i*_M+j)*_L+k];}
    const T & operator()( const int & i, const int & j, const int & k) const {return _data[(i*_M+j)*_L+k];}

    int size() const {return _size;}
    int size1() const {return _N;}
    int size2() const {return _M;}
    int size3() const {return _L;}

    const Vec3d<T> & operator=(const T & a)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=a;
        }
        return *this;
    }
    const Vec3d<T> & operator=(const Vec3d<T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec3d<T> & operator=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec3d<T> & operator+=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]+=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec3d<T> & operator-=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]-=u(i);
        }
        return *this;
    }

    const Vec3d<T> & operator*=(const T & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]*=u;
        }
        return *this;
    }

    const Vec3d<T> & operator/=(const T & u)
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
    int _N;
    int _M;
    int _L;
};

template<typename E1, typename E2, typename T>
class Vec3dPlus : public Vec3dExpression<Vec3dPlus<E1,E2,T> ,T>
{
public:
    Vec3dPlus(const E1 & u, const E2 & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)+_v(i);}
    T operator()(const int & i, const int & j, const int & k) const {return _u(i,j,k)+_v((i*_u.size2()+j)*_u.size3()+k);}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
    int size3() const {return _u.size3();}
private:
    const E1 & _u;
    const E2 & _v;
};

template<typename E1, typename E2, typename T>
Vec3dPlus<E1,E2,T> operator+(const Vec3dExpression<E1,T> & u, const Vec3dExpression<E2,T> & v)
{
    return Vec3dPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
Vec3dPlus<E1,E2,T> operator+(const Vec3dExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    return Vec3dPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
Vec3dPlus<E1,E2,T> operator+(const VecExpression<E2,T> & v, const Vec3dExpression<E1,T> & u)
{
    return Vec3dPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
class Vec3dMinus : public Vec3dExpression<Vec3dMinus<E1,E2,T> ,T>
{
public:
    Vec3dMinus(const E1 & u, const E2 & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)-_v(i);}
    T operator()(const int & i, const int & j, const int & k) const {return _u(i,j,k)-_v((i*_u.size2()+j)*_u.size3()+k);}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
    int size3() const {return _u.size3();}
private:
    const E1 & _u;
    const E2 & _v;
};

template<typename E1, typename E2, typename T>
Vec3dMinus<E1,E2,T> operator-(const Vec3dExpression<E1,T> & u, const Vec3dExpression<E2,T> & v)
{
    return Vec3dMinus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
Vec3dMinus<E1,E2,T> operator-(const Vec3dExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    return Vec3dMinus<E1,E2,T>(u,v);
}

template<typename E, typename T>
class Vec3dPlusScalar : public Vec3dExpression<Vec3dPlusScalar<E,T> ,T>
{
public:
    Vec3dPlusScalar(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)+_v;}
    T operator()(const int & i, const int & j, const int & k ) const {return _u(i,j,k)+_v;}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
    int size3() const {return _u.size3();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
Vec3dPlusScalar<E,T> operator+(const Vec3dExpression<E,T> & u, const T & v)
{
    return Vec3dPlusScalar<E,T>(u,v);
}

template<typename E, typename T>
Vec3dPlusScalar<E,T> operator+(const T & u, const Vec3dExpression<E,T> & v)
{
    return Vec3dPlusScalar<E,T>(v,u);
}

template<typename E, typename T>
class Vec3dProductScalar : public Vec3dExpression<Vec3dProductScalar<E,T> ,T>
{
public:
    Vec3dProductScalar(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)*_v;}
    T operator()(const int & i, const int & j, const int & k ) const {return _u(i,j,k)*_v;}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
    int size3() const {return _u.size3();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
Vec3dProductScalar<E,T> operator*(const Vec3dExpression<E,T> & u, const T & v)
{
    return Vec3dProductScalar<E,T>(u,v);
}

template<typename E, typename T>
Vec3dProductScalar<E,T> operator*(const T & u, const Vec3dExpression<E,T> & v)
{
    return Vec3dProductScalar<E,T>(v,u);
}

template<typename E, typename T>
ostream & operator<<(ostream & cout, const Vec3dExpression<E,T> & u)
{
    for(int i=0;i<u.size1();i++)
    {
        for(int j=0;j<u.size2();j++)
        {
            for(int k=0;k<u.size3();k++)
            {
                cout<<u(i,j,k)<<"  ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
    return cout;
}
#endif // VECTOR3d_H
