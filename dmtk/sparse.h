#ifndef __DMTK_SPARSE_H__
#define __DMTK_SPARSE_H__

#include <valarray>
#include <vector>
#include "conj.h"
#include "vector.h"
#include "matrix.h"
#include "sparse_iter.h"

namespace dmtk
{

#include "meta.h"

template<class T>
class Sparse
{
  private:
    size_t num_rows, num_cols;

    std::vector<int> v_index;
    typename std::vector<T> v_data;

  public:
    Sparse(): num_rows(0), num_cols(0){}

    Sparse(size_t minsize):v_index(minsize), v_data(minsize), num_rows(0), num_cols(0){}

    Sparse(size_t minsize, size_t nrows, size_t ncols):v_index(minsize), v_data(minsize), num_rows(nrows), num_cols(ncols){}

    Sparse(size_t nrows, size_t ncols):v_index(1), v_data(1), num_rows(nrows), num_cols(ncols){}


    explicit Sparse(const Sparse<T>& );
    explicit Sparse(const Matrix<T>& m){ store(m); }
    explicit Sparse(const Matrix<T>&, const Matrix<T>&); // tensor product
    explicit Sparse(const Sparse<T>&, const Sparse<T>&); // tensor product

//  Ranges

    sparse_iter<T> operator()(Range, Range);
    csparse_iter<T> operator()(Range, Range) const;

//  Operators

    T operator()(size_t col, size_t row) const;

    Sparse& operator=(const Matrix<T> &m) { store(m); return *this; }

    template<class Expr> 
    Sparse& operator=(const IterExpr<T,Expr> &m) { store(m); return *this; }

//  Methods

    Sparse& store(const Matrix<T>& m);

    template<class Expr> 
    Sparse& store(const IterExpr<T,Expr>& m);

    Sparse& clear();
    float sparsity() { return (float)size()/(float)(num_rows*num_cols); }

    Sparse& redimension(size_t nrows, size_t ncols) { num_rows = nrows; num_cols = ncols; }

    size_t rows() const { return num_rows; }
    size_t cols() const { return num_cols; }

//  Iterators

    typename std::vector<T>::iterator col_data(size_t row);
    std::vector<int>::iterator col_index(size_t row);
    typename std::vector<T>::const_iterator col_data(size_t row) const;
    std::vector<int>::const_iterator col_index(size_t row) const;

    T diagonal(size_t row) const;
    size_t non_zero(size_t row) const;

//  Memory management

    size_t size() const { return v_index.size(); }
    size_t capacity() const { return v_index.capacity(); }
    Sparse& resize(size_t nrows, size_t ncols) 
      { num_rows = nrows; num_cols = ncols; return *this; }
    Sparse& reserve(size_t size) { v_index.reserve(size); v_data.reserve(size); }
    Sparse& reserve(size_t nrows, size_t ncols) { v_index.reserve(nrows*ncols); v_data.reserve(nrows*ncols); };

//  Reference

    ConstRef<T, Sparse<T> > ref() const { return ConstRef<T, Sparse<T> >(*this); }

};

#include "sparse_implement.h"

template<class T>
inline  
Sparse<T>::Sparse(const Sparse<T> &m) 
{
  v_index = m.v_index;
  v_data = m.v_data;
  num_rows = m.num_rows;
  num_cols = m.num_cols;
}

template<class T>
inline Sparse<T>& 
Sparse<T>::clear() 
{
  v_index.clear(); 
  v_data.clear(); 
  num_rows = num_cols = 0;
  return *this;
}

template<class T>
inline T 
Sparse<T>::operator()(size_t col, size_t row) const 
{
  std::vector<int>::const_iterator irow, icol;
  typename std::vector<T>::const_iterator data;
  int non_zero;
  int pos = 0;

  irow = v_index.begin();
  data = v_data.begin();

  for(pos = 0; pos < row; pos++){
    int step = *irow + 1;
    irow += step; 
    data += step; 
  }
  if(row == col) return *data;

  non_zero = *irow;
  icol = irow + 1;
  data++;

  for(pos = 1; pos <= non_zero; pos++){
    if(*icol == col) return *data;
    icol++; 
    data++; 
  }

  return (T)0;
}

template<class T>
inline typename std::vector<T>::iterator
Sparse<T>::col_data(size_t row) 
{
  std::vector<int>::iterator irow;
  typename std::vector<T>::iterator data;
  int pos = 0, step;

  irow = v_index.begin();
  data = v_data.begin();

  step = *irow++;
  data++;

  for(pos = 0; pos < row; pos++){
    irow += step; 
    data += step; 
    step = *irow++;
    data++;
  }

  return data;
}

template<class T>
inline std::vector<int>::iterator
Sparse<T>::col_index(size_t row) 
{
  std::vector<int>::iterator irow;
  int pos = 0;
  int step;

  irow = v_index.begin();
  step = *irow++;

  for(pos = 0; pos < row; pos++){
    irow += step; 
    step = *irow++;
  }

  return irow;
}

template<class T>
inline typename std::vector<T>::const_iterator
Sparse<T>::col_data(size_t row) const 
{
  std::vector<int>::const_iterator irow;
  typename std::vector<T>::const_iterator data;
  int pos = 0, step;

  irow = v_index.begin();
  data = v_data.begin();

  step = *irow++;
  data++;

  for(pos = 0; pos < row; pos++){
    irow += step; 
    data += step; 
    step = *irow++;
    data++;
  }

  return data;
}

template<class T>
inline std::vector<int>::const_iterator
Sparse<T>::col_index(size_t row) const
{
  std::vector<int>::const_iterator irow;
  int pos = 0;
  int step;

  irow = v_index.begin();
  step = *irow++;

  for(pos = 0; pos < row; pos++){
    irow += step; 
    step = *irow++;
  }

  return irow;
}


template<class T>
inline T
Sparse<T>::diagonal(size_t row) const 
{
  return *(v_data(row)--);
}

template<class T>
inline size_t
Sparse<T>::non_zero(size_t row) const 
{
  return *(v_index(row)--);
}

template<class T>
template<class Expr>
inline Sparse<T>& 
Sparse<T>::store(const IterExpr<T, Expr>& m) 
{
  int non_zero;
  T data;
  int pivot;

  v_index.clear();
  v_data.clear();

  for(int row = 0; row < rows(); row++){
    non_zero = 0;
    pivot = v_index.size();
    v_index.push_back(0);
    v_data.push_back((T)0);

//  Here we add row

    for(int col = 0; col < cols(); col++){
      if((data = m(col,row)) != T(0)){
        if(row == col){
          v_data[pivot] = data;
        }else{
          v_data.push_back(data);
          v_index.push_back(col);
          non_zero++;
        } 
      }
    }
    v_index[pivot] = non_zero;
  }

  return *this;
}

template<class T>
inline Sparse<T>& 
Sparse<T>::store(const Matrix<T> &m) 
{
  int non_zero;
  T data;
  int pivot;

  v_index.clear();
  v_data.clear();

  num_rows = m.rows();
  num_cols = m.cols();

  for(int row = 0; row < m.rows(); row++){
    non_zero = 0;
    pivot = v_index.size();
    v_index.push_back(0);
    v_data.push_back((T)0);

//  Here we add row

    for(int col = 0; col < m.cols(); col++){
      if((data = m[row][col]) != T(0)){
        if(row == col){
          v_data[pivot] = data;
        }else{
          v_data.push_back(data);
          v_index.push_back(col);
          non_zero++;
        } 
      }
    }
    v_index[pivot] = non_zero;
  }

  return *this;
}

template<class T>
inline  
Sparse<T>::Sparse(const Matrix<T>& m1, const Matrix<T>& m2) 
{
  int non_zero;
  T data;
  int row, col;
  int pivot;

  v_index.clear();
  v_data.clear();

  num_rows = m1.rows() * m2.rows();
  num_cols = m1.cols() * m2.cols();

  for(int row1 = 0; row1 < m1.rows(); row1++)
    for(int row2 = 0; row2 < m2.rows(); row2++){
      row = row1 * m2.rows() + row2;
      non_zero = 0;
      pivot = v_index.size();
      v_index.push_back(0);
      v_data.push_back((T)0);

//  Here we add row

      for(int col1 = 0; col1 < m1.cols(); col1++)
        for(int col2 = 0; col2 < m2.cols(); col2++){
          col = col1 * m2.cols() + col2;

          if((data = m1[row1][col1]*m2[row2][col2]) != (T) 0.){
            if(row == col){
              v_data[pivot] = data;
            }else{
              v_data.push_back(data);
              v_index.push_back(col);
              non_zero++;
            } 
          }
        }
        v_index[pivot] = non_zero;
    }

}

template<class T>
inline  
Sparse<T>::Sparse(const Sparse<T>& m1, const Sparse<T>& m2) 
{
  int non_zero;
  T data;
  int row, col;
  int col1, col2;
  int pivot;
  std::vector<int>::const_iterator irow1, irow2, icol1, icol2;
  typename std::vector<T>::const_iterator drow1, drow2, data1, data2;

  v_index.clear();
  v_data.clear();

  num_rows = m1.rows() * m2.rows();
  num_cols = m1.cols() * m2.cols();

  irow1 = m1.v_index.begin();
  drow1 = m1.v_data.begin();
  int nz1 = *irow1;

  for(int row1 = 0; row1 < m1.rows(); row1++){
    irow2 = m2.v_index.begin();
    drow2 = m2.v_data.begin();
    int nz2 = *irow2;

    for(int row2 = 0; row2 < m2.rows(); row2++){

      row = row1 * m2.rows() + row2;

      non_zero = 0;
      pivot = v_index.size();
      v_index.push_back(0);
      v_data.push_back((T)0);

//  Here we add row

      icol1 = irow1;
      data1 = drow1;

      for(int pos1 = 0; pos1 <= nz1; pos1++){

        col1 = pos1 == 0 ? row1 : *icol1; 

        icol2 = irow2;
        data2 = drow2;

        for(int pos2 = 0; pos2 <= nz2; pos2++){

          col2 = pos2 == 0 ? row2 : *icol2; 

          col = col1 * m2.cols() + col2;

	  data = *data1 * *data2;

          if(row == col){
            v_data[pivot] = data;
          }else{
            v_data.push_back(data);
            v_index.push_back(col);
            non_zero++;
          } 
             
          icol2++; 
          data2++; 
        }
        icol1++; 
        data1++; 
      }

      v_index[pivot] = non_zero;

      irow2 += (nz2 + 1);
      drow2 += (nz2 + 1);
    }
    irow1 += (nz1 + 1);
    drow1 += (nz1 + 1);
  }

}

} // namespace dmtk 

#endif // __DMTK_SPARSE_H__ 
