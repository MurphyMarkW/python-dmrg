#ifndef __DMTK_VECTOR_H__
#define __DMTK_VECTOR_H__

#include <iostream>
#include <valarray>
#include <iosfwd>
#include "conj.h"
#include "slice_iter.h"
#include "gslice_iter.h"
#include "range.h"
#include "meta.h"
#include "array_util.h"

namespace dmtk
{

template<class T>
class Vector:public std::valarray<T>
{
  public:
    typedef typename std::valarray<T> valarray;
    typedef T*iterator;
    typedef const T* const_iterator;

    Vector(): 
      valarray(1), v_size(0), v_capacity(1) { init(); };

    explicit Vector(size_t s):
      valarray(s),v_size(s),v_capacity(s) { init(); }

    Vector(size_t s, const T& v):
      valarray(v, s),v_size(s),v_capacity(s) { init(); }

    Vector(size_t s, const T* v):
      valarray(v, s),v_size(s),v_capacity(s) { init(); }

    Vector(const Vector<T>& v):
      valarray(v.capacity()), v_size(v.size()), v_capacity(v.capacity()) 
      { 
        init();
        iterator dest = begin();
        const_iterator orig = v.begin();
        int n = size();
        array_copy(n, orig, dest); 
//        while(n--) *dest++ = *orig++;
      }

    Vector(cslice_iter<T> v):
      valarray(v.size()),v_size(v.size()),v_capacity(v.size()) 
      { 
        init();
        iterator pos = begin();
        int n = std::min(size(), v.size());
        resize(n);
        array_copy2(n, v, *this); 
//        while(n--) *pos++ = *v++;
      }

    Vector(slice_iter<T> v):
      valarray(v.size()),v_size(v.size()),v_capacity(v.size()) 
      { 
        init();
        iterator pos = begin();
        int n = std::min(size(), v.size());
        resize(n);
        array_copy2(n, v, *this); 
//        while(n--) *pos++ = *v++;
      }

    template<class Expr>
    Vector(const IterExpr<T,Expr>& v):
      valarray(v.size2()),v_size(v.size2()),v_capacity(v.size2()) 
      { 
        for(int i = 0; i < size(); i++) this->operator[](i) = v[i];
        init();
      }

//  Operators 

    T operator()(size_t i) const; 
    T& operator()(size_t i);
    T operator[](size_t i) const; // () and [] are equivalent; 
    T& operator[](size_t i);

    Vector& operator=(const Vector<T>&);
    Vector& operator=(slice_iter<T>);
    Vector& operator=(cslice_iter<T>);
    Vector& operator=(const_iterator);
    Vector& operator=(const T&);

    template<class Expr>
    Vector& operator=(const IterExpr<T,Expr>&);

    Vector& operator+=(const Vector<T> &);
    Vector& operator-=(const Vector<T> &);
    Vector& operator*=(const Vector<T> &);
    Vector& operator/=(const Vector<T> &);

    Vector& operator+=(const T &);
    Vector& operator-=(const T &);
    Vector& operator*=(const T &);
    Vector& operator/=(const T &);

    template<class Expr>
    Vector& operator+=(const IterExpr<T,Expr>&);
    template<class Expr>
    Vector& operator-=(const IterExpr<T,Expr>&);
    template<class Expr>
    Vector& operator*=(const IterExpr<T,Expr>&);
    template<class Expr>
    Vector& operator/=(const IterExpr<T,Expr>&);

//  Iterator

    T* array() { return v_start; }
    const T* array() const { return v_start; }
    iterator begin() { return v_start; }
    const_iterator begin() const { return v_start; }
    iterator end() { return v_end; }
    const_iterator end() const { return v_end; }

//  Reference
                                                                                
    ConstRef<T, Vector<T> > ref() const { return ConstRef<T, Vector<T> >(*this); }

//  Ranges

    slice_iter<T> operator()(Range);
    cslice_iter<T> operator()(Range) const;
    gslice_iter<T> as_matrix(size_t ncols,size_t nrows);
    cgslice_iter<T> as_matrix(size_t ncols,size_t nrows) const;

//  Memory management

    inline size_t size() const { return v_size; }
    inline size_t capacity() const { return v_capacity; }
    inline size_t memory() const { return v_capacity * sizeof(T); }
    inline size_t usage() const { return v_size * sizeof(T); }
    Vector& resize(size_t);
    Vector& resize(size_t, const T&);
    Vector& reserve(size_t);

//  IterExpr auxiliary methods
    size_t size1() const { return 1; }
    size_t size2() const { return size(); }

//  Methods

    Vector<T>& randomize(bool positive = false)
    {
      static long idum = abs((long)this)/0xffff;
//      idum = 1234567;
      if(positive){
        for(int i = 0; i < size(); i++)
          operator()(i) = T(quickran(idum));
      } else {
        for(int i = 0; i < size(); i++)
          operator()(i) = T(2.*quickran(idum)-1.);
      }
//      this->operator=(T(1.));
      return *this;
    }

//  Streams

    void read(std::istream& s);
    void write(std::ostream& s) const;

  private:
    T* v_start;
    T* v_end;
    size_t v_size;
    size_t v_capacity; 

    void init() 
      { 
        v_start = &valarray::operator[](0); 
        v_end = &valarray::operator[](v_size); 
      }
};

#include "vector_implement.h"

template<class T>
inline T
Vector<T>::operator()(size_t i) const 
{
  return valarray::operator[](i);
}

template<class T>
inline T&
Vector<T>::operator()(size_t i) 
{
  return valarray::operator[](i);
};

template<class T>
inline T
Vector<T>::operator[](size_t i) const 
{
  return valarray::operator[](i);
}

template<class T>
inline T&
Vector<T>::operator[](size_t i) 
{
  return valarray::operator[](i);
};


template<class T>
slice_iter<T>
Vector<T>::operator()(Range range)
{
  std::slice s(range.start(), range.size(), range.stride());
  slice_iter<T> iter(this, s);
  return iter;
}


template<class T>
cslice_iter<T>
Vector<T>::operator()(Range range) const
{
  std::slice s(range.start(), range.size(), range.stride());
  cslice_iter<T> iter(this, s);
  return iter;
}


template<class T>
inline Vector<T>&
Vector<T>::operator=(const Vector<T> &v)
{
  if(v.size() != size()) resize(v.size());
  int n = size();
  iterator dest = begin();
  const_iterator orig = v.begin();
  array_copy(n, orig, dest);
//  while(n--) *dest++ = *orig++;
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(const T& v)
{
  iterator iter = begin();
  int n = size();
  while(n--) *iter++ = v;
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(const_iterator v)
{
  iterator pos = begin();
  int n = size();
  array_copy(n, v, pos);
//  while(n--) *pos++ = *v++;
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(slice_iter<T> v)
{
  iterator pos = begin();
  int n = v.size();
  resize(n);
  array_copy2(n, v, *this);
//  while(n--) *pos++ = *v++;
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(cslice_iter<T> v)
{
  iterator pos = begin();
  int n = v.size();
  resize(n);
  array_copy2(n, v, *this);
//  while(n--) *pos++ = *v++;
  return *this;
}

template<class T>
template<class Expr>
inline Vector<T>&
Vector<T>::operator=(const IterExpr<T,Expr>& vv)
{
  resize(vv.size2());
  iterator pos = begin();
  for(size_t i = 0; i < size(); i++) *pos++ = vv[i];
//  array_copy2(size(), vv, *this);

  return *this;
}

////////////////////////////////////////////////////////

#define BINARY_OP(op,ap) \
template<class T> \
template<class Expr> \
inline Vector<T>& \
Vector<T>::op(const IterExpr<T,Expr>& v) \
{ \
  iterator pos = begin(); \
  for(int i = 0; pos != end(); i++, pos++) *pos ap v[i]; \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline Vector<T>& \
Vector<T>::op(const Vector<T>& mm) \
{ \
  iterator dest = begin(); \
  const_iterator orig = mm.begin(); \
  int n = size(); \
  while(n--) *dest++ ap *orig++; \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline Vector<T>& \
Vector<T>::op(const T& v) \
{ \
  iterator pos = begin(); \
  while(pos != end()) *pos++ ap v; \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

////////////////////////////////////////////////////////
template<class T>
inline T
product(const Vector<T> &a, const Vector<T> &b)
{
  using namespace std;

  if(a.size() != b.size())
     cerr << "** Warning: Vector sizes do not comform\n";

  typename Vector<T>::const_iterator pa = a.begin();
  typename Vector<T>::const_iterator pb = b.begin();
  int n = std::min(a.size(), b.size());
  return dot_product(n, pa, 1, pb, 1);
}

template<class T>
inline T
dot_product(const Vector<T> &a, const Vector<T> &b)
{
  return product(a,b);
}

template<class T>
Vector<T>&
Vector<T>::resize(size_t n)
{
  if(n > capacity()){ 
     Vector<T> aux(*this);
     T* paux = &aux[0];

     this->valarray::resize(n);
     T* p = &operator[](0);

     array_copy(aux.size(), paux, p);
//     while(paux != aux.end()) *p++ = *paux++;
          
     v_capacity = n;
  }
  v_size = n;   

  init();
  return *this;
}

template<class T>
Vector<T>&
Vector<T>::resize(size_t n, const T& v)
{
  if(n > capacity()){ 
     Vector<T> aux(*this);
     T* paux = &aux[0];

     valarray::resize(n, v);
     T* p = &operator[](0);

     array_copy(aux.size(), paux, p);
//     while(paux != aux.end()) *p++ = *paux++;
          
     v_capacity = n;
  }
  v_size = n;   

  init();
  return *this;
}

template<class T>
Vector<T>&
Vector<T>::reserve(size_t n)
{
  if(n > capacity()){ 
     Vector<T> aux(*this);
     T* paux = &aux[0];

     valarray::resize(n);
     T* p = &operator[](0);

     array_copy(aux.size(), paux, p);
//     while(paux != aux.end()) *p++ = *paux++;
          
     v_capacity = n;
  }

  init();
  return *this;
}

// Streams

template<class T>
std::ostream& operator << (std::ostream& s, const Vector<T>& v)
{
  typename Vector<T>::const_iterator iter = v.begin();
  s << "(" << (*iter++);
  for(; iter != v.end(); iter++)
    s << "," << (*iter);
  s << ")";

  return s; 
}


template<class T>
void 
Vector<T>::write(std::ostream &s) const
{
  size_t n = size();
  s.write((const char *)&n, sizeof(size_t));

  const_iterator iter;
  for(iter = begin(); iter != end(); iter++)
    s.write((const char *)iter, sizeof(T));
}

template<class T>
void 
Vector<T>::read(std::istream &s) 
{
  size_t n;
  s.read((char *)&n, sizeof(size_t));
  resize(n);

  iterator iter;
  for(iter = begin(); iter != end(); iter++)
    s.read((char *)iter, sizeof(T));
}


/////////////////////////////////////////////////////////////////

} // namespace dmtk

#endif // __DMTK_VECTOR_H__ 
