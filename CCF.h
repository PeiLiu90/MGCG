#ifndef CCF_H
#define CCF_H

#include <iostream>
#include <fstream>
using namespace std;

/*matrix coordinates is from 1 to N*/
template<typename T>
class CCF
{
public:
    CCF(){};
    void SetSize(const int & N)
    {
        _size=N;
        _colptr = new int [N+1];
    }
    void SetLength(const int & M)
    {
        _length= M;
        _rowind = new int [M];
        _data = new T [M];
    }
    void SetData(const int * cp, const int * ri, const T * dt)
    {
        for(int i=0;i<=_size;i++)
        {
            _colptr[i] = cp[i];
        }
        for(int i=0;i<_length;i++)
        {
            _rowind[i] = ri[i];
            _data[i] = dt[i];
        }
    }
    void PrintMatrix()
    {
        cout.precision(2);
        fixed(cout);
        for(int i=0;i<_size;i++)
        {
            for(int k= 1;k<_rowind[_colptr[i]-1];k++)
            {
                cout<<"   ";
            }
            cout<<_data[_colptr[i]-1];
            for(int j=_colptr[i];j<_colptr[i+1]-1;j++)
            {
                for(int k=_rowind[j-1]+1;k<_rowind[j];k++)
                {
                    cout<<"   ";
                }
                cout<<"  "<<_data[j];
            }
            for(int k=_rowind[_colptr[i+1]-2];k<_size;k++)
            {
                cout<<"   ";
            }
            cout<<endl;
        }
        ofstream mat("mat.txt");
        mat.precision(2);
        fixed(mat);
        for(int i=0;i<_size;i++)
        {
            for(int k= 1;k<_rowind[_colptr[i]-1];k++)
            {
                mat<<"   ";
            }
            mat<<_data[_colptr[i]-1];
            for(int j=_colptr[i];j<_colptr[i+1]-1;j++)
            {
                for(int k=_rowind[j-1]+1;k<_rowind[j];k++)
                {
                    mat<<"   ";
                }
                mat<<"  "<<_data[j];
            }
            for(int k=_rowind[_colptr[i+1]-2];k<_size;k++)
            {
                mat<<"   ";
            }
            mat<<endl;
        }
    }
    int & colptr(const int & i) {return _colptr[i];}
    const int & colptr(const int & i) const {return _colptr[i];}
    int & rowind(const int &i) {return _rowind[i];}
    const int & rowind(const int & i) const {return _rowind[i];}
    T & data(const int & i) {return _data[i];}
    const T & data(const int & i) const {return _data[i];}
    const int & size() const {return _size;}
private:
    int * _colptr;
    int * _rowind;
    T * _data;
    int _size;
    int _length;
};


#endif // CCF_H
