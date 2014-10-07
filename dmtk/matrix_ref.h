#ifndef __DMTK_MATRIX_REF_H__
#define __DMTK_MATRIX_REF_H__

#include "count_ref.h"

namespace dmtk
{

template<class T>
class MatrixRef:public CountedRef<Matrix<T> >
{
  public:
    MatrixRef() {}
    MatrixRef(const Matrix<T> *p): CountedRef<Matrix<T> >(p) {}

    MatrixRef(const MatrixRef<T>& other) { CountedRef<Matrix<T> >(other); }

    MatrixRef& operator=(const MatrixRef<T>& other) 
      { CountedRef<Matrix<T> >operator=(other); return *this; }

    MatrixRef& operator=(cgslice_iter<T> s);
      { _counted->_pt->operator=(s); return *this; }
    MatrixRef& operator=(gslice_iter<T> s);
      { _counted->_pt->operator=(s); return *this; }
    MatrixRef& operator=(const T *v);
      { _counted->_pt->operator=(v); return *this; }
    MatrixRef& operator=(const T& v);
      { _counted->_pt->operator=(v); return *this; }

    MatrixRef& operator=(const Matrix<T>& other) 
      { _counted->_pt->operator=(other); return *this; }

    template<class Expr>
    MatrixRef& operator=(const IterExpr<T,Expr>& expr)
      { _counted->_pt->operator=(expr); return *this; }

    MatrixRef& operator+=(const MatrixRef<T>& other) 
      { _counted->_pt->operator+=(other); return *this; }
    MatrixRef& operator-=(const MatrixRef<T>& other) 
      { _counted->_pt->operator-=(other); return *this; }
    MatrixRef& operator*=(const MatrixRef<T>& other) 
      { _counted->_pt->operator*=(other); return *this; }
    MatrixRef& operator/=(const MatrixRef<T>& other) 
      { _counted->_pt->operator/=(other); return *this; }

    template<class Expr>
    MatrixRef& operator+=(const IterExpr<T,Expr>& expr)
      { _counted->_pt->operator+=(expr); return *this; }
    template<class Expr>
    MatrixRef& operator-=(const IterExpr<T,Expr>& expr)
      { _counted->_pt->operator-=(expr); return *this; }
    template<class Expr>
    MatrixRef& operator*=(const IterExpr<T,Expr>& expr)
      { _counted->_pt->operator*=(expr); return *this; }
    template<class Expr>
    MatrixRef& operator/=(const IterExpr<T,Expr>& expr)
      { _counted->_pt->operator/=(expr); return *this; }

    const A& operator[](size_t i) const
      { return _counted->_pt->operator[](i); }

    A& operator[](size_t i)
      { return _counted->_pt->operator[](i); }

    const A& operator()(size_t i) const
      { return _counted->_pt->operator()(i); }

    A& operator()(size_t i)
      { return _counted->_pt->operator()(i); }

    const A& operator()(size_t i, size_t j) const
      { return _counted->_pt->operator()(i,j); }

    A& operator()(size_t i, size_t j)
      { return _counted->_pt->operator()(i,j); }

//  Iterators

    slice_iter<T> row(size_t i)  // row iterator
      { return _counted->_pt->row(i); }
    cslice_iter<T> row(size_t i) const  // constant row iterator
      { return _counted->_pt->row(i); }
    slice_iter<T> column(size_t i)  // column iterator
      { return _counted->_pt->column(i); }
    cslice_iter<T> column(size_t i) const  // constant column iterator
      { return _counted->_pt->column(i); }

//  Reference

    ConstRef<T, MatrixRef<T> > ref() const { return ConstRef<T, MatrixRef<T> >(*this);

//  Ranges

    slice_iter<T> operator()(Range r, size_t i);
      { return _counted->_pt->operator()(r,i); }
    cslice_iter<T> operator()(Range r, size_t i) const;
      { return _counted->_pt->operator()(r,i); }
    slice_iter<T> operator()(size_t i, Range r);
      { return _counted->_pt->operator()(i,r); }
    cslice_iter<T> operator()(size_t i, Range r) const;
      { return _counted->_pt->operator()(i,r); }
    gslice_iter<T> operator()(Range r, Range r);
      { return _counted->_pt->operator()(r,r); }
    cgslice_iter<T> operator()(Range r, Range r) const;
      { return _counted->_pt->operator()(r,r); }
    slice_iter<T> diagonal();
      { return _counted->_pt->diagonal(); }
    cslice_iter<T> diagonal() const;
      { return _counted->_pt->diagonal(); }

//  Methods

    inline size_t rows() const { return _counted->_pt->rows(); }
    inline size_t cols() const { return _counted->_pt->cols(); }

//  Memory management

    inline size_t memory() const 
      { return _counted->_pt->diagonal(); }
    inline size_t usage() const 
      { return _counted->_pt->diagonal(); }

    MatrixRef<T>& resize(size_t ncols, size_t nrows)
      { return _counted->_pt->resize(ncols, nrows); }
    MatrixRef<T>& resize(size_t ncols, size_t nrows, const T &v)
      { return _counted->_pt->resize(ncols, nrows, v); }
    MatrixRef<T>& reserve(size_t ncols, size_t nrows)
      { return _counted->_pt->reserve(ncols, nrows); }
 
}




} // namespace dmtk

#endif // __DMTK_MATRIX_REF_H__
