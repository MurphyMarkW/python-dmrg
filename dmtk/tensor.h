#ifndef __DMTK_STATE_H__
#define __DMTK_STATE_H__

#include <iosfwd>
#include "vector.h"

using namespace std;

namespace dmtk
{

#include "meta.h"

template <class T>
class Tensor: public Vector<T>
{
  public:

    typedef std::vector<StateSpace> _V;
    typedef std::vector<StateSpace>::iterator iterator;
    typedef std::vector<StateSpace>::const_iterator const_iterator;

    Tensor(): {}
    Tensor(size_t _d): Vector<T>(_d){}
    Tensor(size_t _d1, size_t _d2, size_t _d3, size_t _d4): Vector<T>(_d1*_d2*_d3*_d4) {}
    Tensor(const Vector<T> &v) Vector<T>(v) {}
    Tensor(const Tensor<T>& v): Vector<T>(v) {}

    Tensor& operator=(const T& v)
      { Vector<T>::operator=(v); return *this; }
    
    Tensor& operator=(const Tensor<T>& v)
      { Vector<T>::operator=(v); return *this; }

    Tensor& operator=(const Vector<T>& v)
      { Vector<T>::operator=(v); return *this; };

    template<class Expr>
    Tensor& operator=(const IterExpr<T,Expr>& exp)
      { Vector<T>::operator=(exp); return *this; }

    Tensor& operator+=(const Tensor<T> &v)
      { Vector<T>::operator+=(v); return *this; }
    Tensor& operator-=(const Tensor<T> &v)
      { Vector<T>::operator-=(v); return *this; }
    Tensor& operator*=(const Tensor<T> &v)
      { Vector<T>::operator*=(v); return *this; }
    Tensor& operator/=(const Tensor<T> &v)
      { Vector<T>::operator/=(v); return *this; }

    Tensor& operator+=(const T &v)
      { Vector<T>::operator+=(v); return *this; }
    Tensor& operator-=(const T &v)
      { Vector<T>::operator-=(v); return *this; }
    Tensor& operator*=(const T &v)
      { Vector<T>::operator*=(v); return *this; }
    Tensor& operator/=(const T &v)
      { Vector<T>::operator/=(v); return *this; }

    template<class Expr>
    Tensor& operator+=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator+=(exp); return *this; }
    template<class Expr>
    Tensor& operator-=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator-=(exp); return *this; }
    template<class Expr>
    Tensor& operator*=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator*=(exp); return *this; }
    template<class Expr>
    Tensor& operator/=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator/=(exp); return *this; }

    T& operator[](size_t i)
      { return Vector<T>::operator[](i); }

    T operator[](size_t i) const
      { return Vector<T>::operator[](i); }

    T& operator()(size_t i)
      { return Vector<T>::operator[](i); }

    T operator()(size_t i) const
      { return Vector<T>::operator[](i); }


    T& operator()(size_t i1, size_t i2, size_t i3, size_t i4)
      {
        return Vector<T>::operator[](i1*d2*d3*d4+i2*d3*d4+i3*d4+i4);
      }

    T operator()(size_t i1, size_t i2, size_t i3, size_t i4) const
      {
        return Vector<T>::operator[](i1*d2*d3*d4+i2*d3*d4+i3*d4+i4);
      }

    T& operator()(size_t i1, size_t i2, size_t i3, size_t i4, int mask)
      {
        if(mask & MASK_BLOCK1 && mask & MASK_BLOCK2)
          return operator()(i1,i2,i3,i4);
        else if(mask & MASK_BLOCK2 && mask & MASK_BLOCK3)
          return operator()(i3,i1,i2,i4);
        else if(mask & MASK_BLOCK3 && mask & MASK_BLOCK4)
          return operator()(i3,i4,i1,i2);
      }

    Tensor& resize(size_t d1, size_t d2, size_t d3, size_t d4) 
      { 
         Vector<T>::resize(d1*d*d3*d4);
         return *this;
      } 

//  Streams

    void write(std::ostream &s) const
    {
      Vector<T>::write(s);
    }

    void read(std::istream &s) 
    {
      Vector<T>::read(s);
    }

};


////////////////////////////////////////////////////////

} // namespace dmtk

#endif // __DMTK_STATE_H__
