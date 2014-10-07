#ifndef __DMTK_SPARSE_IMPLEMENT_H__
#define __DMTK_SPARSE_IMPLEMENT_H__

/////////////////////////////////////////////////////////////////

template<class T>
sparse_iter<T>
Sparse<T>::operator()(Range col_range, Range row_range)
{
  slice s_col(col_range.start(),
           col_range.size(),
           col_range.stride());
  slice s_row(row_range.start(),
           row_range.size(),
           row_range.stride());
  sparse_iter<T> iter(this, s_col, s_row);
  return iter;
}

template<class T>
csparse_iter<T>
Sparse<T>::operator()(Range col_range, Range row_range) const
{
  slice s_col(col_range.start(),
           col_range.size(),
           col_range.stride());
  slice s_row(row_range.start(),
           row_range.size(),
           row_range.stride());
  csparse_iter<T> iter(this, s_col, s_row);
  return iter;
}

/////////////////////////////////////////////////////////////////

#define SPARSE_SPARSE_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, ConstRef<T,Sparse<T> >, ConstRef<T,Sparse<T> >, ap<T, ConstRef<T,Sparse<T> >, ConstRef<T,Sparse<T> > > > > \
op(Sparse<T>& a, Sparse<T>& b) \
{ \
  typedef ConstRef<T,Sparse<T> > IT; \
  typedef IterBinOp<T, IT, IT, ap<T, IT, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b.ref())); \
}

SPARSE_SPARSE_OP(operator+, DMApAdd)
SPARSE_SPARSE_OP(operator-, DMApSubs)
SPARSE_SPARSE_OP(operator*, DMApMul)
SPARSE_SPARSE_OP(operator/, DMApDiv)
#undef SPARSE_SPARSE_OP

/////////////////////////////////////////////////////////////////

#define EXPR_SPARSE_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, ConstRef<T,Sparse<T> >, ap<T, IterExpr<T, Expr>, ConstRef<T,Sparse<T> > > > > \
op(IterExpr<T, Expr> a, Sparse<T>& b) \
{ \
  typedef ConstRef<T,Sparse<T> > IT; \
  typedef IterBinOp<T, IterExpr<T, Expr>, IT, ap<T, IterExpr<T, Expr>, IT > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.ref())); \
}

EXPR_SPARSE_OP(operator+, DMApAdd)
EXPR_SPARSE_OP(operator-, DMApSubs)
EXPR_SPARSE_OP(operator*, DMApMul)
EXPR_SPARSE_OP(operator/, DMApDiv)
#undef EXPR_SPARSE_OP

/////////////////////////////////////////////////////////////////

#define SPARSE_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, ConstRef<T,Sparse<T> >, IterExpr<T, Expr> , ap<T, ConstRef<T,Sparse<T> >, IterExpr<T, Expr> > > > \
op(Sparse<T>& a, IterExpr<T, Expr> b) \
{ \
  typedef ConstRef<T,Sparse<T> > IT; \
  typedef IterBinOp<T, IT, IterExpr<T, Expr>, ap<T, IT, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b)); \
}

SPARSE_EXPR_OP(operator+, DMApAdd)
SPARSE_EXPR_OP(operator-, DMApSubs)
SPARSE_EXPR_OP(operator*, DMApMul)
SPARSE_EXPR_OP(operator/, DMApDiv)
#undef SPARSE_EXPR_OP

/////////////////////////////////////////////////////////////////

#define SPARSE_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, ConstRef<T,Sparse<T> >, T, ap<T, ConstRef<T,Sparse<T> >, T> > > \
op(Sparse<T>& a, const T& s) \
{ \
  typedef ConstRef<T,Sparse<T> > IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), s)); \
}

SPARSE_SCALAR_OP(operator+, DMApAdd)
SPARSE_SCALAR_OP(operator-, DMApSubs)
SPARSE_SCALAR_OP(operator*, DMApMul)
SPARSE_SCALAR_OP(operator/, DMApDiv)
#undef SPARSE_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_SPARSE_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, ConstRef<T,Sparse<T> >, T, ap<T, ConstRef<T,Sparse<T> >, T> > > \
op(const T& s, Sparse<T>& a) \
{ \
  typedef ConstRef<T,Sparse<T> > IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), s)); \
}

SCALAR_SPARSE_OP(operator+, DMApAdd)
SCALAR_SPARSE_OP(operator-, DMApSubs)
SCALAR_SPARSE_OP(operator*, DMApMul)
SCALAR_SPARSE_OP(operator/, DMApDiv)
#undef SCALAR_SPARSE_OP

/////////////////////////////////////////////////////////////////

template<class T>
class SparseSparseTensor
{
  public:
    static inline T apply(ConstRef<T,Sparse<T> > m1, ConstRef<T,Sparse<T> > m2, size_t i, size_t j)
    {
      const Sparse<T> &s1 = m1.ref();
      const Sparse<T> &s2 = m2.ref();
      size_t col1 = i / s2.cols();
      size_t col2 = i - col1 * s2.cols();
      size_t row1 = j / s2.rows();
      size_t row2 = j - row1 * s2.rows();
      return m1(col1,row1)*m2(col2,row2);
    }
};

template<class T>
IterExpr<T, IterBinOp<T, ConstRef<T,Sparse<T> >, ConstRef<T,Sparse<T> >, SparseSparseTensor<T> > >
tensor(const Sparse<T>& a, const Sparse<T>& b)
{
  typedef ConstRef<T,Sparse<T> > IT;
  typedef IterBinOp<T, IT, IT, SparseSparseTensor<T> > BinExpr;
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b.ref()));
}

/////////////////////////////////////////////////////////////////

template<class T>
class SparseTranspose
{
  public:
    static inline T apply(ConstRef<T,Sparse<T> > m, size_t i, size_t j)
    {
      return m(j, i);
    }
};

template<class T>
class SparseConj
{
  public:
    static inline T apply(ConstRef<T,Sparse<T> > m, size_t i, size_t j)
    {
      return std::conj(m(i, j));
    }
};

#define SPARSE_METHOD(op, ap) \
template<class T> \
IterExpr<T, IterOp<T, ConstRef<T,Sparse<T> >, ap<T> > > \
op(const Sparse<T> &a) \
{ \
  typedef ConstRef<T,Sparse<T> > IT; \
  typedef IterOp<T, IT, ap<T> > Expr; \
  return IterExpr<T, Expr > (Expr(a.ref())); \
}

SPARSE_METHOD(transpose, SparseTranspose);
SPARSE_METHOD(conj, SparseConj);
#undef SPARSE_METHOD

/////////////////////////////////////////////////////////////////

template<class T>
class SparseVectorProduct
{
  public:
    static inline T apply(ConstRef<T,Sparse<T> > m, typename Vector<T>::const_iterator v, size_t n)
    {
      const Sparse<T> &s = m.ref();
      std::vector<int>::const_iterator irow, icol;
      typename std::vector<T>::const_iterator data;
      int non_zero;
      int pos;
      T res = T(0);

      irow = s.col_index(0);
      data = s.col_data(0);
      irow--; // set to origin
      data--;

      for(pos = 0; pos < n; pos++){
        int step = *irow + 1;
        irow += step;
        data += step;
      }

      non_zero = *irow++;
      res += *data++ * v[n];

      icol = irow;
      for(pos = 1; pos <= non_zero; pos++){
        int col = *icol++;
        res += *data++ * v[col];
      }

      return res;
    }
};

template<class T>
IterExpr<T, IterBinOp<T, ConstRef<T,Sparse<T> >, typename Vector<T>::const_iterator, SparseVectorProduct<T> > >
product(const Sparse<T>& m, const Vector<T>& v)
{
  typedef IterBinOp<T, ConstRef<T,Sparse<T> >, typename Vector<T>::const_iterator, SparseVectorProduct<T> > BinExpr;

  return IterExpr<T, BinExpr > (BinExpr(m.ref(), v.begin()));
}

/////////////////////////////////////////////////////////////////

template<class T>
class SparseSparseProduct
{
  public:
    static inline T apply(ConstRef<T,Sparse<T> > m1, ConstRef<T,Sparse<T> > m2, size_t i, size_t j)
    {
      const Sparse<T> &s1 = m1.ref();
      const Sparse<T> &s2 = m2.ref();
      T res = T(0);

      for(size_t pos = 0; pos < s1.cols(); pos++)
        res += s1(j, pos) * s2(pos, i);
    }
};

template<class T>
IterExpr<T, IterBinOp<T, ConstRef<T,Sparse<T> >, ConstRef<T,Sparse<T> >, SparseSparseProduct<T> > >
product(const Sparse<T>& m1, const Sparse<T>& m2)
{
  typedef ConstRef<T,Sparse<T> > IT;
  typedef IterBinOp<T, IT, IT, SparseSparseProduct<T> > BinExpr;
  return IterExpr<T, BinExpr > (BinExpr(m1.ref(), m2.ref()));
}

/////////////////////////////////////////////////////////////////

#endif /* __DMTK_SPARSE_IMPLEMENT_H__ */
