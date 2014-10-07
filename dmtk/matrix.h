#ifndef __DMTK_MATRIX_H__
#define __DMTK_MATRIX_H__

// Generic Matrix class using templates
// Matrix is a valarray subclass
// The data is stored in Fortran style: column0.column1....column(n-1) 
// for compatibility with numeric software like BLAS and LAPACK
// To retrieve the data in Fortran style use m(col,row)
// To retrieve the data in C style use m[row][col]

#include <valarray>
#include <vector>
#include <iosfwd>
#include <iostream>
#include "conj.h"
#include "slice_iter.h"
#include "gslice_iter.h"
#include "vector.h"
#include "range.h"
#include "array_util.h"

namespace dmtk
{

#include "meta.h"

template<class T>class Sparse; // defined later in sparse.h

template<class T>
class Matrix:public std::valarray<T>
{
  private:
    T* _data;
    size_t num_rows, num_cols;
    size_t v_capacity;

    void init() { _data = &valarray::operator[](0); }

    Matrix& matrix_from_sparse(const Sparse<T>&);
  public:
    typedef typename std::valarray<T> valarray;
    typedef T value_type;

//    Matrix():valarray(4),num_rows(2),num_cols(2),v_capacity(4) { init(); }
    Matrix(): valarray(1), num_rows(1), num_cols(1), v_capacity(1) { init(); }

    Matrix(size_t ncols, size_t nrows):valarray(nrows*ncols),num_rows(nrows),num_cols(ncols),v_capacity(nrows * ncols) { init(); }
    Matrix(size_t ncols, size_t nrows, const T& v):valarray(v, nrows*ncols),num_rows(nrows),num_cols(ncols),v_capacity(nrows * ncols) { init(); }

    Matrix(size_t ncols, size_t nrows, const T* v):valarray(v, nrows*ncols),num_rows(nrows),num_cols(ncols),v_capacity(nrows * ncols) { init(); }

    Matrix(const Matrix<T>& m):valarray(m.rows()*m.cols()),num_rows(m.rows()),num_cols(m.cols()),v_capacity(m.v_capacity)
      { 
        init(); 
        const T* orig = m.array();
        T *dest = array();

        size_t n = num_rows*num_cols;
        array_copy(n, orig, dest);
//        while(n--) *dest++ = *orig++; 
      }

    Matrix(const Sparse<T>& s): valarray(s.cols()*s.rows()),num_rows(s.rows()),num_cols(s.cols()),v_capacity(s.rows()*s.cols()) { matrix_from_sparse(s); }

    Matrix(Vector<T>&, Vector<T>&); // tensor product of two vectors
				    // c_ij = a_i * b_j

    Matrix(Matrix<T>&, Matrix<T>&); // tensor product of two matrices 
				    // c_ijkl = a_ij * b_kl

    Matrix(cgslice_iter<T> s): valarray(s.size1()*s.size2()),num_rows(s.size2()),num_cols(s.size1()),v_capacity(s.size1()*s.size2())
      {
        for(size_t i = 0; i < rows(); i++)
          for(size_t j = 0; j < cols(); j++)
            this->operator()(j,i) = s(j,i);

        init();
      }

    Matrix(gslice_iter<T> s): valarray(s.size1()*s.size2()),num_rows(s.size2()),num_cols(s.size1()),v_capacity(s.size1()*s.size2())
      {
        for(size_t i = 0; i < rows(); i++)
          for(size_t j = 0; j < cols(); j++)
            this->operator()(j,i) = s(j,i);

        init();
      }

    template<class Expr>
    Matrix(const IterExpr<T,Expr>& s):valarray(s.size1()*s.size2()),num_rows(s.size2()),num_cols(s.size1()),v_capacity(s.size1()*s.size2())
      {
        for(size_t i = 0; i < rows(); i++)
          for(size_t j = 0; j < cols(); j++)
            this->operator()(j,i) = s(j,i);

        init();
      }
//  Operators

    T operator()(size_t col, size_t row) const; // Fortran style
    T& operator()(size_t col, size_t row);
    slice_iter<T> operator[] (size_t row); // C style
    cslice_iter<T> operator[] (size_t row) const; // C style

//  Asignment 
 
    Matrix& operator=(const Matrix<T> &);
    Matrix& operator=(const Sparse<T>& s) 
      { matrix_from_sparse(s); return *this; }

    Matrix& operator=(cgslice_iter<T> s);
    Matrix& operator=(gslice_iter<T> s);
    Matrix& operator=(const T *v);
    Matrix& operator=(const T& v);
    template<class Expr>
    Matrix& operator=(const IterExpr<T,Expr>&);

    Matrix& operator+=(const Matrix<T> &);
    Matrix& operator-=(const Matrix<T> &);
    Matrix& operator*=(const Matrix<T> &);
    Matrix& operator/=(const Matrix<T> &);

    Matrix& operator+=(const T &);
    Matrix& operator-=(const T &);
    Matrix& operator*=(const T &);
    Matrix& operator/=(const T &);

    template<class Expr>
    Matrix& operator+=(const IterExpr<T,Expr>&);
    template<class Expr>
    Matrix& operator-=(const IterExpr<T,Expr>&);
    template<class Expr>
    Matrix& operator*=(const IterExpr<T,Expr>&);
    template<class Expr>
    Matrix& operator/=(const IterExpr<T,Expr>&);

//  Iterators

    slice_iter<T> row(size_t);  // row iterator
    cslice_iter<T> row(size_t) const ;  // constant row iterator
    slice_iter<T> column(size_t);  // column iterator
    cslice_iter<T> column(size_t) const;  // constant column iterator
    T* array() { return _data; }
    const T* array() const { return _data; }

//  Reference

    ConstRef<T, Matrix<T> > ref() const { return ConstRef<T, Matrix<T> >(*this); }

//  Ranges

    slice_iter<T> operator()(Range, size_t);
    cslice_iter<T> operator()(Range, size_t) const;
    slice_iter<T> operator()(size_t, Range);
    cslice_iter<T> operator()(size_t, Range) const;
    gslice_iter<T> operator()(Range, Range);
    cgslice_iter<T> operator()(Range, Range) const;
    slice_iter<T> diagonal();
    cslice_iter<T> diagonal() const;
    slice_iter<T> as_vector();
    cslice_iter<T> as_vector() const;

//  Methods
   
    Matrix<T>& rotate(const Matrix<T>& u);
    Matrix<T> t() const;
    Matrix<T> ct() const;
    Matrix<T>& diagonalize(Vector<double> &w);
    Matrix<T>& diagonalize(Vector<float> &w);
    Matrix<T>& randomize(int seed = -1)
    {
      static long idum = seed == -1 ? abs((long)this)/0xffff : seed;
      for(int i = 0; i < cols(); i++)
        for(int j = 0; j < rows(); j++)
          operator()(i,j) = quickran(idum);
      return *this;
    }

    Matrix<T>& invert();
    T det() const;
    inline T trace() const 
    { 
      T r = T(0);
      for(uint i = 0; i < std::min(num_rows, num_cols); i++) 
        r += operator()(i,i);
      return r;
    }

    inline size_t rows() const { return num_rows; }
    inline size_t cols() const { return num_cols; }

//  Memory management

    inline size_t memory() const { return v_capacity * sizeof(T); }
    inline size_t usage() const { return num_rows * num_cols * sizeof(T); }
    Matrix<T>& reshape(size_t ncols, size_t nrows);
    Matrix<T>& resize(size_t ncols, size_t nrows);
    Matrix<T>& resize(size_t ncols, size_t nrows, const T &);
    Matrix<T>& reserve(size_t ncols, size_t nrows);

//  IterExpr auxiliary methods
    inline size_t size1() const { return num_cols; }
    inline size_t size2() const { return num_rows; }

//  Streams

    void read(std::istream& s);
    void write(std::ostream& s) const;
};

#include "matrix_implement.h" 

//////////////////////////////////////////////////////////

// Contructor

template<class T>
inline 
Matrix<T>::Matrix(Vector<T> &v1, Vector<T>& v2)
{
  valarray(v1.size() * v2.size());
  *this = tensor(v1,v2);
  init();
}

template<class T>
inline 
Matrix<T>::Matrix(Matrix<T> &m1, Matrix<T>& m2)
{
  valarray(m1.rows()*m2.rows()*m1.cols()*m2.cols());
  *this = tensor(m1, m2);
  init();
}

// Assignment

template<class T>
inline Matrix<T>& 
Matrix<T>::operator=(const Matrix<T> &m) 
{
  if(rows() != m.num_rows || cols() != m.num_cols)
    reshape(m.num_cols, m.num_rows);

  const T* orig = m.array();
  T *dest = array();

  size_t n = num_rows*num_cols;
  array_copy(n, orig, dest);
//  while(n--) *dest++ = *orig++; 

  return *this;
}

template<class T>
inline Matrix<T>& 
Matrix<T>::operator=(const T *v) 
{
  T *dest = array();
  size_t n = num_rows*num_cols;
  array_copy(n, v, dest);
//  while(n--) *dest++ = *v++; 

  return *this;
}

template<class T>
inline Matrix<T>& 
Matrix<T>::operator=(const T& v) 
{
  T *dest = array();
  size_t n = num_rows*num_cols;
  while(n--) *dest++ = v; 
  return *this;
}

template<class T>
inline Matrix<T>&
Matrix<T>::operator=(cgslice_iter<T> mm)
{
  reshape(mm.size1(), mm.size2());
  for(size_t i = 0; i < rows(); i++)
    for(size_t j = 0; j < cols(); j++)
      this->operator()(j,i) = mm(j,i);

  return *this;
}

template<class T>
inline Matrix<T>&
Matrix<T>::operator=(gslice_iter<T> mm)
{
  reshape(mm.size1(), mm.size2());
  for(size_t i = 0; i < rows(); i++)
    for(size_t j = 0; j < cols(); j++)
      this->operator()(j,i) = mm(j,i);

  return *this;
}

template<class T>
template<class Expr>
inline Matrix<T>&
Matrix<T>::operator=(const IterExpr<T,Expr>& mm)
{
  if(mm.size1() != 0 && mm.size2() != 0) reshape(mm.size1(), mm.size2());
  for(size_t i = 0; i < rows(); i++)
    for(size_t j = 0; j < cols(); j++){
      operator()(j,i) = mm(j,i);
    }

  return *this;
}

//////////////////////////////////////////////////////////

#define BINARY_OP(op,ap) \
template<class T> \
template<class Expr> \
inline Matrix<T>& \
Matrix<T>::op(const IterExpr<T,Expr>& mm) \
{ \
  for(size_t i = 0; i < rows(); i++) \
    for(size_t j = 0; j < cols(); j++) \
      operator()(j,i) ap mm(j,i); \
 \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline Matrix<T>& \
Matrix<T>::op(const Matrix<T>& mm) \
{ \
  reshape(mm.cols(),mm.rows()); \
  const T* orig = mm.array(); \
  T* dest = array(); \
  int n = mm.cols()*mm.rows(); \
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
inline Matrix<T>& \
Matrix<T>::op(const T& v) \
{ \
  for(size_t i = 0; i < rows(); i++) \
    for(size_t j = 0; j < cols(); j++) \
      operator()(j,i) ap v; \
 \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op) \
//////////////////////////////////////////////////////////

// Reference

template<class T>
inline slice_iter<T> 
Matrix<T>::row(size_t i) 
{
  return slice_iter<T>(this, std::slice(i, num_cols, num_rows));
}

template<class T>
inline cslice_iter<T> 
Matrix<T>::row(size_t i) const 
{
  return cslice_iter<T>(this, std::slice(i, num_cols, num_rows));
}

template<class T>
inline slice_iter<T> 
Matrix<T>::column(size_t i) 
{
  return slice_iter<T>(this, std::slice(i*num_rows, num_rows, 1));
}

template<class T>
inline cslice_iter<T> 
Matrix<T>::column(size_t i) const
{
  return cslice_iter<T>(this, std::slice(i*num_rows, num_rows, 1));
}

template<class T>
inline T 
Matrix<T>::operator()(size_t ncol, size_t nrow) const 
{
  return valarray::operator[](ncol*num_rows+nrow);
}

template<class T>
inline T&
Matrix<T>::operator()(size_t ncol, size_t nrow) 
{
  return valarray::operator[](ncol*num_rows+nrow);
}

template<class T>
inline slice_iter<T> 
Matrix<T>::operator[](size_t nrow) 
{
  return row(nrow);
}

template<class T>
inline cslice_iter<T> 
Matrix<T>::operator[](size_t nrow) const 
{
  return row(nrow);
}

//////////////////////////////////////////////////////////

// Memory management

template<class T>
Matrix<T>&
Matrix<T>::resize(size_t ncols, size_t nrows)
{
  if(nrows == num_rows && ncols == num_cols) return *this;

  Matrix<T> aux(*this);
  size_t new_size = nrows * ncols;
  T *orig, *dest;

  if(new_size > v_capacity){ 
    valarray::resize(new_size);
    v_capacity = new_size;
  } 

  for(uint col = 0 ; col < std::min(num_cols, ncols); col++){
    orig = &aux.column(col)[0];
    dest = &column(col)[0];
    for(uint row = 0 ; row < std::min(num_rows, nrows) ; row++){
        *dest++ = *orig++; 
    }
  }
  num_cols = ncols;
  num_rows = nrows;

  _data = &valarray::operator[](0);

  return *this; 
}

template<class T>
Matrix<T>&
Matrix<T>::resize(size_t ncols, size_t nrows, const T& v)
{
  if(nrows == num_rows && ncols == num_cols) return *this;

  Matrix<T> aux(*this);
  T *orig, *dest;
  size_t old_size = num_rows * num_cols;
  size_t new_size = nrows * ncols;

  if(new_size > v_capacity){ 
    valarray::resize(new_size, v);
    v_capacity = new_size;
  }

  for(uint col = 0 ; col < std::min(num_cols, ncols); col++){
    orig = &aux.column(col)[0];
    dest = &column(col)[0];
    for(uint row = 0 ; row < std::min(num_rows, nrows) ; row++){
        *dest++ = *orig++; 
    }
  }
  num_cols = ncols;
  num_rows = nrows;

  _data = &valarray::operator[](0);

  return *this; 
}

template<class T>
Matrix<T>&
Matrix<T>::reshape(size_t ncols, size_t nrows)
{
  if(nrows == num_rows && ncols == num_cols) return *this;

  size_t new_size = nrows * ncols;

  if(new_size > v_capacity){ 
    valarray::resize(new_size);
    v_capacity = new_size;
  }

  num_cols = ncols;
  num_rows = nrows;

  _data = &valarray::operator[](0);

  return *this; 
}

template<class T>
Matrix<T>&
Matrix<T>::reserve(size_t ncols, size_t nrows)
{
  Matrix<T> aux(*this);
  T *orig, *dest;
  size_t pos;
  size_t old_size = num_rows * num_cols;
  size_t new_size = nrows * ncols;

  if(new_size > v_capacity){ 
    valarray::resize(new_size);
    v_capacity = new_size;

    for(uint col = 0 ; col < num_cols; col++){
      orig = &aux.column(col)[0];
      dest = &column(col)[0];
      for(uint row = 0 ; row < num_rows ; row++){
          *dest++ = *orig++; 
      }
    }
  }

  _data = &valarray::operator[](0);

  return *this; 
}


//////////////////////////////////////////////////////////

template<class T>
inline Matrix<T>
Matrix<T>::t() const
{
  Matrix<T> aux(cols(),rows());

  aux = transpose(*this);
  return aux; 
}

template<class T>
inline Matrix<T>
Matrix<T>::ct() const
{
  Matrix<T> aux(cols(),rows());

  aux = ctranspose(*this);
  return aux; 
}

template<class T>
inline Matrix<T>&
Matrix<T>::rotate(const Matrix& u)
{
  Matrix<T> aux(u), aux2(u);

  aux = product(*this,u);
  aux2 = ctranspose(u);
  *this = product(aux2, aux);
 
  return *this; 
}

template<class T>
inline Matrix<T>&
Matrix<T>::invert()
{
  using namespace std;

  if(num_rows != num_cols)
    cerr << "** Warning: Matrix dimensions are different\n";

  Matrix<T> aux(*this);
  int l, j, i, m = std::min(num_rows, num_cols);
  T r;

  for(i = 0; i < m; i++){
    r = T(1) / aux(i,i);
    aux(i,i) = T(0);
    for(j = 0; j < m; j++){
      aux(j,i) *= r;
      for(l = 0; l < m; l++) aux(j,l) = aux(j,l) - aux(i,l)*aux(j,i);
    }
    for(j = 0; j < m; j++) aux(i,j) *= (-r);
    aux(i,i) = r;
  }

  *this = aux;
  return *this;
}

template<class T>
inline T
Matrix<T>::det() const
{
  using namespace std;

  if(num_rows != num_cols)
    cerr << "** Warning: Matrix dimensions are different\n";

  Matrix<T> *aux(*this);
  int k, l, j, i, m = std::min(num_rows, num_cols);
  T r;

  for(i = 0; i < m; i++){
    r = T(1) / aux(i,i);
    aux(i,i) = T(0);
    for(j = 0; j < m; j++){
      aux(j,i) *= r;
      for(l = 0; l < m; l++) aux(j,l) = aux(j,l) - aux(i,l)*aux(j,i);
    }
    for(j = 0; j < m; j++) aux(i,j) *= (-r);
    aux(i,i) = r;
  }

  T d = T(1);
  for(k = 0; k < m; k++) d *= aux(k,k);

  return d;
}

template<class T>
inline Matrix<T>&
Matrix<T>::matrix_from_sparse(const Sparse<T>& s)
{
  using namespace std;

  vector<int>::const_iterator irow, icol;
  typename vector<T>::const_iterator data;
  int col, pos;

  using namespace std;

  if(num_cols != s.cols() || num_rows != s.rows()){
    cerr << "** Warning: Matrix and Sparse sizes do not comform\n";
    cerr << num_cols << " " << num_rows << " " << s.cols() << " " << s.rows() << endl;
  }

  resize(s.cols(), s.rows());

  irow = s.col_index(0);
  data = s.col_data(0);
  irow--; // set to origin
  data--;

  for(uint row = 0; row < rows(); row++){
    int non_zero = *irow;
    icol = irow++;
    icol++;

    this->operator()(row, row) = *data++;
    for(int pos = 1; pos <= non_zero; pos++){
      col = *icol++;
      this->operator()(col, row) = *data++;
    }
    irow = icol;
  }
  return *this;
}

////////////////////////////////////////////////////////
template<>
Matrix<float>& 
Matrix<float>::diagonalize(Vector<float>& ev)
{
  Matrix<float> &b = *this;
  ev.resize(b.rows()); // eigenvalues

#ifdef WITH_LAPACK
  Vector<float> d(b.rows());
  Vector<float> work(3*b.rows());
  DMTK_int info = 0;

  DMTK_int b_rows = b.rows();
  DMTK_int w_size = work.size();
  ssyev_('V','U',b_rows,b.array(),b_rows, 
         d.array(),work.array(),w_size,info);

  if(info != 0) cout << "*** ERROR in dsyev: info != 0 (failed to converge)\n";

  ev = d;

#else // !WITH_LAPACK
//***** NOTE: This is to use when Lapack is not available

  Vector<float> e(b.rows()+1);
  Vector<float> d(b.rows()+1);
  Matrix<float> z(b.rows()+1, b.rows()+1);

  z(Range(1,b.rows()),Range(1,b.rows())) = b(Range(0,b.rows()-1),Range(0,b.rows()-1));

  tred2(z, b.rows(), d, e, true);
  tqli(d, e, b.rows(), z, true);

  b(Range(0,b.rows()-1),Range(0,b.rows()-1)) = z(Range(1,b.rows()),Range(1,b.rows()));
  ev(Range(0,b.rows()-1)) = d(Range(1,b.rows()));

#endif // WITH LAPACK

  return *this;
}

template<>
Matrix<double>& 
Matrix<double>::diagonalize(Vector<double>& ev)
{
  Matrix<double> &b = *this;
  ev.resize(b.rows()); // eigenvalues

#ifdef WITH_LAPACK
  Vector<double> d(b.rows());
  Vector<double> work(3*b.rows());
  DMTK_int info = 0;

  DMTK_int b_rows = b.rows();
  DMTK_int w_size = work.size();
  dsyev_('V','U',b_rows,b.array(),b_rows, 
         d.array(),work.array(),w_size,info);

  if(info != 0) cout << "*** ERROR in dsyev: info != 0 (failed to converge)\n";

  ev = d;

#else // !WITH_LAPACK
//***** NOTE: This is to use when Lapack is not available

  Vector<double> e(b.rows()+1);
  Vector<double> d(b.rows()+1);
  Matrix<double> z(b.rows()+1, b.rows()+1);

  z(Range(1,b.rows()),Range(1,b.rows())) = b(Range(0,b.rows()-1),Range(0,b.rows()-1));

  tred2(z, b.rows(), d, e, true);
  tqli(d, e, b.rows(), z, true);

  b(Range(0,b.rows()-1),Range(0,b.rows()-1)) = z(Range(1,b.rows()),Range(1,b.rows()));
  ev(Range(0,b.rows()-1)) = d(Range(1,b.rows()));

#endif // WITH LAPACK

  return *this;
}

template<>
Matrix<complex<double> >& 
Matrix<complex<double> >::diagonalize(Vector<double>& ev)
{
  Matrix<complex<double> > &b = *this;
  ev.resize(b.rows()); // eigenvalues

#ifdef WITH_LAPACK

  Vector<double> d(b.rows());
  Vector<complex<double> > work(2*b.rows());
  Vector<double> rwork(3*b.rows()-2);
  DMTK_int info = 0;
  DMTK_int b_rows = b.rows();
  DMTK_int w_size = work.size();
  zheev_('V','U',b_rows,b.array(),b_rows, 
         d.array(),work.array(),w_size,rwork.array(),info);

  if(info != 0) cout << "*** ERROR in zheev: info != 0 (failed to converge)\n";

  ev = d;

#else // !WITH_LAPACK
//------ NOTE: This part uses TQLI with a matrix twice the size of the original
//       (see Numerical Recipes). However, there is a bug somewhere. Use Lapack.

  Vector<double> e(2*b.rows()+1);
  Vector<double> d(2*b.rows()+1);
  Matrix<double> z(2*b.rows()+1, 2*b.rows()+1);

  for(int i = 1; i <= b.rows(); i++)
    for(int j = 1; j <= b.rows(); j++){
      z(i,j) = real(b(i-1,j-1));
      z(b.cols()+i,b.rows()+j) = real(b(i-1,j-1));
      z(b.cols()+i,j) = -imag(b(i-1,j-1));
      z(i,b.rows()+i) = imag(b(i-1,j-1));
    }

  tred2(z, 2*b.rows(), d, e, true);
  tqli(d, e, 2*b.rows(), z, true);

//  for(int i = 1; i <= 2*b.rows(); i++) cout << i << " " << d(i) << endl;

// Store eigenstates, checking for orthogonality first

  int n=0;
  for(int i = 1; i <= 2*b.rows(); i++){
    Vector<complex<double> > coef(b.rows());
    for(int j = 1; j <= b.rows(); j++) 
      coef(j-1) = complex<double>(z(i,j),z(i,j+b.rows()));

    bool found = false;
    complex<double> xdeg;
    for(int i1 = 0; i1 < n; i1++){
      xdeg = complex<double>(0.f,0.f);
      for(int i2 = 0; i2 < b.rows(); i2++) xdeg += coef(i2)*conj(b(i1,i2));
      if(real(xdeg*conj(xdeg)) > 1.e-3) { found = true; break; }
//      if(real(xdeg*conj(xdeg)) > 1.e-3) { found = true; nr++; cout << nr << " " << xdeg*conj(xdeg) << endl; break; }
    } 
    if(!found){
      b(n,Range(0,b.rows()-1)) = coef(Range(0,b.rows()-1));
      ev(n) = d(i);
//      xsum += d(i);
      n++;
    }
  }

#endif // WITH_LAPACK
  return *this;
}



////////////////////////////////////////////////////////
template<class T>
inline T
dot_product(const Matrix<T> &a, const Matrix<T> &b)
{
  using namespace std;
                                                                                
  if(a.size() != b.size())
     cerr << "** Warning: Matrix sizes do not comform\n";
                                                                                
  const T* pa = a.array();
  const T* pb = b.array();
                               
  int n = std::min(a.size(), b.size());
  return dot_product(n, pa, 1, pb, 1);
}

template<class T>
Matrix<T>
cexp(const Matrix<T> &m)
{
  if(m.rows() != m.cols()){
     cerr << "** Warning: Matrix should be square\n";
     return m;
  }

// We diagonalize and exponentiate
  Vector<double> wk(m.cols());
  Matrix<T> u = m;
  u.diagonalize(wk);

//  cout << "---------------------------\n";

  int i = 0;
  Matrix<complex<double> > ut(u.rows(),u.cols());
  ut = ctranspose(u);

  Matrix<T> auxm(m);
  auxm = complex<double>(0.,0.);
  for(int j = 0; j < u.rows(); j++){
    auxm(j,j) = std::exp(complex<double>(0.,wk[j]));
  } 
  Matrix<complex<double> > aux(u.cols(),m.rows());
  aux = product(auxm,ut);
  auxm = product(u,aux);

/*
  for(int l = 0; l < m.rows(); l++)
   for(int j = 0; j < m.rows(); j++) cout << l << " " << j << " " << m(l,j) << endl;
*/

  return auxm;
}

template<class T>
Matrix<T>
exp(const Matrix<T> &m)
{
  if(m.rows() != m.cols()){
     cerr << "** Warning: Matrix should be square\n";
     return m;
  }

// We diagonalize and exponentiate
  Vector<double> wk(m.cols());
  Matrix<T> u = m;
  u.diagonalize(wk);

//  cout << "---------------------------\n";

  int i = 0;
//  for(biter = rotated_hij.begin(); biter != rotated_hij.end(); biter++){
  Matrix<complex<double> > ut(u.rows(),u.cols());
  ut = ctranspose(u);

  Matrix<T> auxm(m);
  auxm = complex<double>(0.,0.);
  for(int j = 0; j < u.rows(); j++){
    auxm(j,j) = std::exp(wk[j]);
  } 
  Matrix<complex<double> > aux(u.cols(),m.rows());
  aux = product(auxm,ut);
  auxm = product(u,aux);

/*
  for(int l = 0; l < m.rows(); l++)
   for(int j = 0; j < m.rows(); j++) cout << l << " " << j << " " << m(l,j) << endl;
*/

  return auxm;
}
///////////////////////////////////////////////////////////////////////////

template<class T>
void
Matrix<T>::write(std::ostream &s) const
{
  s.write((const char *)&num_rows, sizeof(size_t));
  s.write((const char *)&num_cols, sizeof(size_t));

  for(uint col = 0; col < cols(); col++) 
    for(uint row = 0; row < rows(); row++){ 
      T v = operator()(col,row);
      s.write((const char *)&v, sizeof(T));
    }
}

template<class T>
void
Matrix<T>::read(std::istream &s)
{
  size_t nr, nc;
  s.read((char *)&nr, sizeof(size_t));
  s.read((char *)&nc, sizeof(size_t));
  resize(nc,nr);

  for(uint col = 0; col < cols(); col++) 
    for(uint row = 0; row < rows(); row++){ 
      T& v = operator()(col,row);
      s.read((char *)&v, sizeof(T));
    }
}

} // namespace dmtk 

#endif // __DMTK_MATRIX_H__ 
