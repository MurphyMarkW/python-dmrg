#ifndef __DMTK_SPARSE_ITER_IMPLEMENT_H__
#define __DMTK_SPARSE_ITER_IMPLEMENT_H__

/////////////////////////////////////////////////////////////////

#define SPARSE_SPARSE_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, csparse_iter<T>, csparse_iter<T>, ap<T, csparse_iter<T>, csparse_iter<T> > > > \
op(csparse_iter<T>& a, csparse_iter<T>& b) \
{ \
  typedef csparse_iter<T> IT; \
  typedef IterBinOp<T, IT, IT, ap<T, IT, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b.begin())); \
}

SPARSE_SPARSE_OP(operator+, DMApAdd)
SPARSE_SPARSE_OP(operator-, DMApSubs)
SPARSE_SPARSE_OP(operator*, DMApMul)
SPARSE_SPARSE_OP(operator/, DMApDiv)
#undef SPARSE_SPARSE_OP

/////////////////////////////////////////////////////////////////

#define EXPR_SPARSE_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, csparse_iter<T>, ap<T, IterExpr<T, Expr>, csparse_iter<T> > > > \
op(IterExpr<T, Expr> a, csparse_iter<T>& b) \
{ \
  typedef csparse_iter<T> IT; \
  typedef IterBinOp<T, IterExpr<T, Expr>, IT, ap<T, IterExpr<T, Expr>, IT > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.begin())); \
}

EXPR_SPARSE_OP(operator+, DMApAdd)
EXPR_SPARSE_OP(operator-, DMApSubs)
EXPR_SPARSE_OP(operator*, DMApMul)
EXPR_SPARSE_OP(operator/, DMApDiv)
#undef EXPR_SPARSE_OP

/////////////////////////////////////////////////////////////////

#define SPARSE_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, csparse_iter<T>, IterExpr<T, Expr> , ap<T, csparse_iter<T>, IterExpr<T, Expr> > > > \
op(csparse_iter<T>& a, IterExpr<T, Expr> b) \
{ \
  typedef csparse_iter<T> IT; \
  typedef IterBinOp<T, IT, IterExpr<T, Expr>, ap<T, IT, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b)); \
}

SPARSE_EXPR_OP(operator+, DMApAdd)
SPARSE_EXPR_OP(operator-, DMApSubs)
SPARSE_EXPR_OP(operator*, DMApMul)
SPARSE_EXPR_OP(operator/, DMApDiv)
#undef SPARSE_EXPR_OP

/////////////////////////////////////////////////////////////////

#define SPARSE_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, csparse_iter<T>, T, ap<T, csparse_iter<T>, T> > > \
op(csparse_iter<T>& a, const T& s) \
{ \
  typedef csparse_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SPARSE_SCALAR_OP(operator+, DMApAdd)
SPARSE_SCALAR_OP(operator-, DMApSubs)
SPARSE_SCALAR_OP(operator*, DMApMul)
SPARSE_SCALAR_OP(operator/, DMApDiv)
#undef SPARSE_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_SPARSE_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, csparse_iter<T>, T, ap<T, csparse_iter<T>, T> > > \
op(const T& s, csparse_iter<T>& a) \
{ \
  typedef csparse_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SCALAR_SPARSE_OP(operator+, DMApAdd)
SCALAR_SPARSE_OP(operator-, DMApSubs)
SCALAR_SPARSE_OP(operator*, DMApMul)
SCALAR_SPARSE_OP(operator/, DMApDiv)
#undef SCALAR_SPARSE_OP

/////////////////////////////////////////////////////////////////

#endif /* __DMTK_SPARSE_ITER_IMPLEMENT_H__ */
