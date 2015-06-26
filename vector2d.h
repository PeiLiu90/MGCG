#ifndef vector2d_H
#define vector2d_H

#include "math.h"
#include "vec.h"

template<typename E, typename T>
class Vec2dExpression : public VecExpression<E,T>
{
public:
    T operator()(const int & i) const{return static_cast<const E&>(*this)(i);}
    T operator()(const int & i , const int & j) const {return static_cast<const E&>(*this)(i,j);}
    int size() const {return static_cast<const E &>(*this).size();}
    int size1() const {return static_cast<const E &>(*this).size1();}
    int size2() const {return static_cast<const E &>(*this).size2();}

    operator E&() { return static_cast<E&>(*this);}
    operator E const&() const { return static_cast<const E&>(*this);}
};

template<typename T>
class Vec2d: public Vec2dExpression<Vec2d<T>,T>, public vec<T>
{
public:
    Vec2d(){};
    Vec2d(const Vec2d<T> & u)
    {
        _N=u.size1();
        _M=u.size2();
        _size= _N*_M;
        _data = new T [_size];
        for(int i=0;i<_N*_M;i++)
        {
            _data[i]=u(i);
        }
    }
    template<typename E>
    Vec2d(const Vec2dExpression<E,T> & u)
    {
        _N= u.size1();
        _M= u.size2();
        _size= _N*_M;
        _data = new T [_size];
        for(int i=0;i<_N*_M;i++)
        {
            _data[i]=u(i);
        }
    }
    Vec2d(const int & n, const int & m) {_N=n;_M=m;_size= _N*_M; _data = new T [_size];}
    void SetSize(const int & N, const int & M){_N= N;_M= M;_size= _N*_M; _data = new T [_size];}

    T & operator()(const int &i) {return _data[i];}
    const T & operator()(const int &i) const {return _data[i];}

    T & operator() (const int & i, const int & j) {return _data[i*_M+j];}
    const T & operator()( const int & i, const int & j) const {return _data[i*_M+j];}

    int size() const {return _N*_M;}
    int size1() const {return _N;}
    int size2() const {return _M;}

    const Vec2d<T> & operator=(const Vec2d<T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec2d<T> & operator=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec2d<T> & operator+=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]+=u(i);
        }
        return *this;
    }

    template<typename E>
    const Vec2d<T> & operator-=(const VecExpression<E,T> & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]-=u(i);
        }
        return *this;
    }

    const Vec2d<T> & operator*=(const T & u)
    {
        for(int i=0;i<_size;i++)
        {
            _data[i]*=u;
        }
        return *this;
    }

    const Vec2d<T> & operator/=(const T & u)
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
};

template<typename E1, typename E2, typename T>
class Vec2dPlus : public Vec2dExpression<Vec2dPlus<E1,E2,T> ,T>
{
public:
    Vec2dPlus(const E1 & u, const E2 & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)+_v(i);}
    T operator()(const int & i, const int & j ) const {return _u(i,j)+_v(i*_u.size2()+j);}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
private:
    const E1 & _u;
    const E2 & _v;
};

template<typename E1, typename E2, typename T>
Vec2dPlus<E1,E2,T> operator+(const Vec2dExpression<E1,T> & u, const Vec2dExpression<E2,T> & v)
{
    return Vec2dPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
Vec2dPlus<E1,E2,T> operator+(const Vec2dExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    return Vec2dPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
Vec2dPlus<E1,E2,T> operator+(const VecExpression<E2,T> & v, const Vec2dExpression<E1,T> & u)
{
    return Vec2dPlus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
class Vec2dMinus : public Vec2dExpression<Vec2dMinus<E1,E2,T> ,T>
{
public:
    Vec2dMinus(const E1 & u, const E2 & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)-_v(i);}
    T operator()(const int & i, const int & j ) const {return _u(i,j)-_v(i*_u.size2()+j);}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
private:
    const E1 & _u;
    const E2 & _v;
};

template<typename E1, typename E2, typename T>
Vec2dMinus<E1,E2,T> operator-(const Vec2dExpression<E1,T> & u, const Vec2dExpression<E2,T> & v)
{
    return Vec2dMinus<E1,E2,T>(u,v);
}

template<typename E1, typename E2, typename T>
Vec2dMinus<E1,E2,T> operator-(const Vec2dExpression<E1,T> & u, const VecExpression<E2,T> & v)
{
    return Vec2dMinus<E1,E2,T>(u,v);
}

template<typename E, typename T>
class Vec2dPlusScalar : public Vec2dExpression<Vec2dPlusScalar<E,T> ,T>
{
public:
    Vec2dPlusScalar(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)+_v;}
    T operator()(const int & i, const int & j ) const {return _u(i,j)+_v;}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
Vec2dPlusScalar<E,T> operator+(const Vec2dExpression<E,T> & u, const T & v)
{
    return Vec2dPlusScalar<E,T>(u,v);
}

template<typename E, typename T>
Vec2dPlusScalar<E,T> operator+(const T & u, const Vec2dExpression<E,T> & v)
{
    return Vec2dPlusScalar<E,T>(v,u);
}

template<typename E, typename T>
class Vec2dProductScalar : public Vec2dExpression<Vec2dProductScalar<E,T> ,T>
{
public:
    Vec2dProductScalar(const E & u, const T & v):_u(u),_v(v){};
    T operator()(const int & i) const { return _u(i)*_v;}
    T operator()(const int & i, const int & j ) const {return _u(i,j)*_v;}
    int size() const {return _u.size();}
    int size1() const {return _u.size1();}
    int size2() const {return _u.size2();}
private:
    const E & _u;
    const T & _v;
};

template<typename E, typename T>
Vec2dProductScalar<E,T> operator*(const Vec2dExpression<E,T> & u, const T & v)
{
    return Vec2dProductScalar<E,T>(u,v);
}

template<typename E, typename T>
Vec2dProductScalar<E,T> operator*(const T & u, const Vec2dExpression<E,T> & v)
{
    return Vec2dProductScalar<E,T>(v,u);
}

template<typename E, typename T>
ostream & operator<<(ostream & cout, const Vec2dExpression<E,T> & u)
{
    for(int i=0;i<u.size1();i++)
    {
        for(int j=0;j<u.size2();j++)
        {
            cout<<u(i,j)<<"  ";
        }
        cout<<endl;
    }
    return cout;
}

#endif // vector2d_H
