#ifndef __DMTK_SPARSE_ITER_H__
#define __DMTK_SPARSE_ITER_H__

// gslice iterator class. From Stroustrup's book.

namespace dmtk
{

#include "meta.h"

template <class T> class Sparse; // defined in sparse.h

template <class T>
class sparse_iter
{
  public:
    sparse_iter(Sparse<T>* vv, slice s1, slice s2):
      v(vv),s1(s1),s2(s2),curr1(0),curr2(0){};
    sparse_iter(const sparse_iter<T>& ss):
      v(ss.v),s1(ss.s1),s2(ss.s2),curr1(ss.curr1),curr2(ss.curr2){}

    sparse_iter begin()  
      {sparse_iter t = *this; curr1 = 0; curr2 = 0; return t; }
    sparse_iter end()  
      {sparse_iter t = *this; curr1 = s1.size(); curr2 = s2.size(); return t; }

    sparse_iter& operator=(const sparse_iter<T>& ss)
      { 
        size_t n1 = std::min(ss.s1.size(), s1.size());
        size_t n2 = std::min(ss.s2.size(), s2.size());
        for(size_t i = 0; i < n1; i++)
          for(size_t j = 0; j < n2; j++)
            this->operator()(i,j) = ss(i,j);
        return *this; 
      }

//  Operators

    size_t size() const { return s1.size(); }
    sparse_iter& operator++() { curr1++; return *this; }
    sparse_iter  operator++(int) {sparse_iter t = *this; curr1++; return t;}

    T& operator()(size_t i, size_t j)
      {return ref(curr1=i, curr2=j);}
    T operator()(size_t i, size_t j) const 
      {return ref(curr1=i, curr2=j);}

  private:
    Sparse<T>* v;
    slice s1;
    slice s2;
    size_t curr1, curr2;

    T& ref(size_t i, size_t j) 
      { 
        return (*v)(s1.start()+i*s1.stride(),s2.start()+j*s2.stride()); 
      }
};

template <class T>
class csparse_iter
{
  public:
    csparse_iter(const Sparse<T>* vv, slice s1, slice s2):
      v(vv),s1(s1),s2(s2),curr1(0),curr2(0){};
    csparse_iter(const csparse_iter<T>& ss):
      v(ss.v),s1(ss.s1),s2(ss.s2),curr1(ss.curr1),curr2(ss.curr2){}

    csparse_iter begin()  
      {csparse_iter t = *this; curr1 = 0; curr2 = 0; return t; }
    csparse_iter end()  
      {csparse_iter t = *this; curr1 = s1.size(); curr2 = s2.size(); return t; }

    csparse_iter& operator=(const csparse_iter<T>& ss)
      { 
        v = ss.v; s1 = ss.s1, s2 = ss.s2; curr1 = ss.curr1; curr2 = ss.curr2;
        return *this; 
      }

//  Operators

    size_t size() const { return s1.size(); }
    csparse_iter& operator++() { curr1++; return *this; }
    csparse_iter operator++(int) {csparse_iter t = *this; curr1++; return t;}

    T operator()(size_t i, size_t j)
      {return ref(curr1=i, curr2=j);}

  private:
    const Sparse<T>* v;
    slice s1;
    slice s2;
    size_t curr1, curr2;

    T ref(size_t i, size_t j) 
      { 
        return (*v)(s1.start()+i*s1.stride(),s2.start()+j*s2.stride()); 
      }
};

#include "sparse_iter_implement.h"

} // namespace dmtk 

#endif // __DMTK_SPARSE_ITER_H__


