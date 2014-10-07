#ifndef __DMTK_LA_UTIL_H__
#define __DMTK_LA_UTIL_H__

#include <complex>
#include "conj.h"
#include "lapack_interface.h"

#ifdef WITH_CBLAS
extern "C" {
#include <atlas_misc.h>
#include <atlas_cblastypealias.h>
#include <atlas_level1.h>
#include <atlas_level2.h>
#include <atlas_level3.h>
#include <cblas.h>
#ifdef WITH_PTHREADS
   #include "atlas_ptalias1.h"
   #include "atlas_ptalias2.h"
   #include "atlas_ptalias3.h"
#endif
}
#endif // WITH_CBLAS

using namespace std;

namespace dmtk
{

#ifdef WITH_LAPACK

void
matrix_matrix_product(const char &c1, const char &c2, 
                      const dmtk::Matrix<float> &a, 
                      const dmtk::Matrix<float> &b, 
                      dmtk::Matrix<float> &res,
                      const float &coefa = 1.,
                      const float &coefb = 0.)
{
  DMTK_int m = c1 == 'N' ? a.rows() : a.cols();
  DMTK_int n = c2 == 'N' ? b.cols() : b.rows();
  DMTK_int k = c1 == 'N' ? a.cols() : a.rows();
  DMTK_int a_rows = a.rows();
  DMTK_int b_rows = b.rows();
  DMTK_int res_rows = res.rows();

#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t1 = CblasNoTrans;
  CBLAS_TRANSPOSE t2 = CblasNoTrans;
  if(c1 == 'C') 
    t1 = CblasConjTrans;
  else if(c1 == 'T')
    t1 = CblasTrans;
  if(c2 == 'C') 
    t2 = CblasConjTrans;
  else if(c2 == 'T')
    t2 = CblasTrans;
//  cblas_dgemm(CblasColMajor, t1, t2, m, n, k, coefa, a.array(), a.rows(), b.array(), b.rows(), coefb, res.array(), res.rows());    
  ATL_sgemm(t1, t2, m, n, k, coefa, a.array(), a_rows, b.array(), b_rows, coefb, res.array(), res_rows);    
#else
  sgemm_(c1, c2, m, n, k, coefa, a.array(), a_rows, b.array(), b_rows, coefb, res.array(), res_rows);    
#endif
}

void
matrix_matrix_product(const char &c1, const char &c2, 
                      const dmtk::Matrix<double> &a, 
                      const dmtk::Matrix<double> &b, 
                      dmtk::Matrix<double> &res,
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  DMTK_int m = c1 == 'N' ? a.rows() : a.cols();
  DMTK_int n = c2 == 'N' ? b.cols() : b.rows();
  DMTK_int k = c1 == 'N' ? a.cols() : a.rows();
  DMTK_int a_rows = a.rows(); 
  DMTK_int b_rows = b.rows(); 
  DMTK_int res_rows = res.rows();

#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t1 = CblasNoTrans;
  CBLAS_TRANSPOSE t2 = CblasNoTrans;
  if(c1 == 'C') 
    t1 = CblasConjTrans;
  else if(c1 == 'T')
    t1 = CblasTrans;
  if(c2 == 'C') 
    t2 = CblasConjTrans;
  else if(c2 == 'T')
    t2 = CblasTrans;
//  cblas_dgemm(CblasColMajor, t1, t2, m, n, k, coefa, a.array(), a.rows(), b.array(), b.rows(), coefb, res.array(), res.rows());    
  ATL_dgemm(t1, t2, m, n, k, coefa, a.array(), a_rows, b.array(), b_rows, coefb, res.array(), res_rows);    
#else
  dgemm_(c1, c2, m, n, k, coefa, a.array(), a_rows, b.array(), b_rows, coefb, res.array(), res_rows);    
#endif
}

void
matrix_matrix_product(const char &c1, const char &c2, 
                      const dmtk::Matrix<complex<double> >& a, 
                      const dmtk::Matrix<complex<double> >& b, 
                      dmtk::Matrix<complex<double> >& res,
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  DMTK_int m = c1 == 'N' ? a.rows() : a.cols();
  DMTK_int n = c2 == 'N' ? b.cols() : b.rows();
  DMTK_int k = c1 == 'N' ? a.cols() : a.rows();
  DMTK_int a_rows = a.rows();
  DMTK_int b_rows = b.rows();
  DMTK_int res_rows = res.rows();

#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t1 = CblasNoTrans;
  CBLAS_TRANSPOSE t2 = CblasNoTrans;
  if(c1 == 'C') 
    t1 = CblasConjTrans;
  else if(c1 == 'T')
    t1 = CblasTrans;
  if(c2 == 'C') 
    t2 = CblasConjTrans;
  else if(c2 == 'T')
    t2 = CblasTrans;
//  cblas_zgemm(CblasColMajor, t1, t2, m, n, k, (const void *)&coefa, (const void *)a.array(), a.rows(), (const void *)b.array(), b.rows(), (const void *)&aux, (void *)res.array(), res.rows());    
  ATL_zgemm(t1, t2, m, n, k, (const double *)&coefa, (const double *)a.array(), a_rows, (const double *)b.array(), b_rows, (const double *)&coefb, (double *)res.array(), res_rows);    
#else
  zgemm_(c1, c2, m, n, k, coefa, a.array(), a_rows, b.array(), b_rows, coefb, res.array(), res_rows);    
#endif
}

void
matrix_matrix_product(const char &c1, const char &c2, 
                      const int &_m,
                      const int &_n,
                      const int &_k,
                      const float *a, 
		      const int &_lda,
                      const float *b, 
		      const int &_ldb,
                      float *res,
		      const int &_ldc,
                      const float &coefa = 1.,
                      const float &coefb = 0.)
{
  DMTK_int m = _m;
  DMTK_int n = _n;
  DMTK_int k = _k;
  DMTK_int lda = _lda;
  DMTK_int ldb = _ldb;
  DMTK_int ldc = _ldc;
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t1 = CblasNoTrans;
  CBLAS_TRANSPOSE t2 = CblasNoTrans;
  if(c1 == 'C') 
    t1 = CblasConjTrans;
  else if(c1 == 'T')
    t1 = CblasTrans;
  if(c2 == 'C') 
    t2 = CblasConjTrans;
  else if(c2 == 'T')
    t2 = CblasTrans;
//  cblas_dgemm(CblasColMajor, t1, t2, m, n, k, coef, a, lda, b, ldb, 0.f, res, ldc);    
  ATL_sgemm(t1, t2, m, n, k, coefa, a, lda, b, ldb, coefb, res, ldc);    
#else
  sgemm_(c1, c2, m, n, k, coefa, a, lda, b, ldb, coefb, res, ldc);    
#endif
}


void
matrix_matrix_product(const char &c1, const char &c2, 
                      const int &_m,
                      const int &_n,
                      const int &_k,
                      const double *a, 
		      const int &_lda,
                      const double *b, 
		      const int &_ldb,
                      double *res,
		      const int &_ldc,
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  DMTK_int m = _m;
  DMTK_int n = _n;
  DMTK_int k = _k;
  DMTK_int lda = _lda;
  DMTK_int ldb = _ldb;
  DMTK_int ldc = _ldc;
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t1 = CblasNoTrans;
  CBLAS_TRANSPOSE t2 = CblasNoTrans;
  if(c1 == 'C') 
    t1 = CblasConjTrans;
  else if(c1 == 'T')
    t1 = CblasTrans;
  if(c2 == 'C') 
    t2 = CblasConjTrans;
  else if(c2 == 'T')
    t2 = CblasTrans;
//  cblas_dgemm(CblasColMajor, t1, t2, m, n, k, coef, a, lda, b, ldb, 0.f, res, ldc);    
  ATL_dgemm(t1, t2, m, n, k, coefa, a, lda, b, ldb, coefb, res, ldc);    
#else
  dgemm_(c1, c2, m, n, k, coefa, a, lda, b, ldb, coefb, res, ldc);    
#endif
}

void
matrix_matrix_product(const char &c1, const char &c2, 
                      const int &_m,
                      const int &_n,
                      const int &_k,
                      const complex<double>* a, 
                      const int &_lda,
                      const complex<double>* b, 
                      const int &_ldb,
                      complex<double>* res,
                      const int &_ldc,
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  DMTK_int m = _m;
  DMTK_int n = _n;
  DMTK_int k = _k;
  DMTK_int lda = _lda;
  DMTK_int ldb = _ldb;
  DMTK_int ldc = _ldc;
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t1 = CblasNoTrans;
  CBLAS_TRANSPOSE t2 = CblasNoTrans;
  if(c1 == 'C') 
    t1 = CblasConjTrans;
  else if(c1 == 'T')
    t1 = CblasTrans;
  if(c2 == 'C') 
    t2 = CblasConjTrans;
  else if(c2 == 'T')
    t2 = CblasTrans;
//  cblas_zgemm(CblasColMajor, t1, t2, m, n, k, (const void *)&coef, (const void *)a, lda, (const void *)b, ldb, (const void *)&coefb, (void *)res, ldc);    
  ATL_zgemm(t1, t2, m, n, k, (const double *)&coefa, (const double *)a, lda, (const double *)b, ldb, (const double *)&coefb, (double *)res, ldc);    
#else
  zgemm_(c1, c2, m, n, k, coefa, a, lda, b, ldb, coefb, res, ldc);    
#endif
}


void
matrix_vector_product(const char& c,
                      const dmtk::Matrix<float> &a, 
                      const dmtk::Vector<float> &v, 
                      dmtk::Vector<float> &res, 
                      const float &coefa = 1.f,
                      const float &coefb = 0.f)
{
  DMTK_int a_rows = a.rows();
  DMTK_int a_cols = a.cols();
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t = CblasNoTrans;
  if(c == 'C') 
    t = CblasConjTrans;
  else if(c == 'T')
    t = CblasTrans;
  cblas_sgemv(CblasColMajor, t, a_rows, a_cols, coefa, a.array(), a_rows, v.array(), 1, coefb, res.array(), 1);    
#else
  sgemv_(c, a_rows, a_cols, coefa, a.array(), a_rows, v.array(), 1, coefb, res.array(), 1);    
#endif
}


void
matrix_vector_product(const char& c,
                      const dmtk::Matrix<double> &a, 
                      const dmtk::Vector<double> &v, 
                      dmtk::Vector<double> &res, 
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  DMTK_int a_rows = a.rows();
  DMTK_int a_cols = a.cols();
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t = CblasNoTrans;
  if(c == 'C') 
    t = CblasConjTrans;
  else if(c == 'T')
    t = CblasTrans;
  cblas_dgemv(CblasColMajor, t, a_rows, a_cols, coefa, a.array(), a_rows, v.array(), 1, coefb, res.array(), 1);    
#else
  dgemv_(c, a_rows, a_cols, coefa, a.array(), a_rows, v.array(), 1, coefb, res.array(), 1);    
#endif
}

void
matrix_vector_product(const char& c,
                      const dmtk::Matrix<complex<double> >&a, 
                      const dmtk::Vector<complex<double> >&v, 
                      dmtk::Vector<complex<double> >&res, 
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  DMTK_int a_rows = a.rows();
  DMTK_int a_cols = a.cols();
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t = CblasNoTrans;
  if(c == 'C') 
    t = CblasConjTrans;
  else if(c == 'T')
    t = CblasTrans;
  cblas_zgemv(CblasColMajor, t, a_rows, a_cols, (const void *)&coefa, a.array(), a_rows, (const void *)v.array(), 1, (const void *)&coefb, (void *)res.array(), 1);    
#else
  zgemv_(c, a_rows, a_cols, coefa, a.array(), a_rows, v.array(), 1, coefb, res.array(), 1);    
#endif
}

void
matrix_vector_product(const char& c,
		      const int& _m,
		      const int& _n,
                      const float *a, 
                      const float *v, 
                      float *res, 
                      const float &coefa = 1.,
                      const float &coefb = 0.)
{
  DMTK_int m = _m;
  DMTK_int n = _n;
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t = CblasNoTrans;
  if(c == 'C') 
    t = CblasConjTrans;
  else if(c == 'T')
    t = CblasTrans;
  cblas_sgemv(CblasColMajor, t, m, n, coefa, a, m, v, 1, coefb, res, 1);    
#else
  sgemv_(c, m, n, coefa, a, m, v, 1, coefb, res, 1);    
#endif
}


void
matrix_vector_product(const char& c,
		      const int& _m,
		      const int& _n,
                      const double *a, 
                      const double *v, 
                      double *res, 
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  DMTK_int m = _m;
  DMTK_int n = _n;
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t = CblasNoTrans;
  if(c == 'C') 
    t = CblasConjTrans;
  else if(c == 'T')
    t = CblasTrans;
  cblas_dgemv(CblasColMajor, t, m, n, coefa, a, m, v, 1, coefb, res, 1);    
#else
  dgemv_(c, m, n, coefa, a, m, v, 1, coefb, res, 1);    
#endif
}

void
matrix_vector_product(const char& c,
		      const int& _m,
		      const int& _n,
                      const complex<double> *a, 
                      const complex<double> *v, 
                      complex<double> *res, 
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  DMTK_int m = _m;
  DMTK_int n = _n;
#ifdef WITH_CBLAS
  CBLAS_TRANSPOSE t = CblasNoTrans;
  if(c == 'C') 
    t = CblasConjTrans;
  else if(c == 'T')
    t = CblasTrans;
  cblas_zgemv(CblasColMajor, t, m, n, (const void *)&coefa, a, m, (const void *)v, 1, (const void *)&coefb, (void *)res, 1);    
#else
  zgemv_(c, m, n, coefa, a, m, v, 1, coefb, res, 1);    
#endif
}

void
vector_scalar_product(Vector<float> &y, const float& a, Vector<float>&x)
{
  DMTK_int size = y.size();
  saxpy_(size, a, x.array(), 1, y.array(), 1);
}

void
vector_scalar_product(Vector<double> &y, const double& a, Vector<double>&x)
{
  DMTK_int size = y.size();
  daxpy_(size, a, x.array(), 1, y.array(), 1);
}

#else // !WITH_LAPACK

template<class T>
void
matrix_matrix_product(const char &c1, const char& c2,
                      const dmtk::Matrix<T>& a, const dmtk::Matrix<T> &b, 
		      dmtk::Matrix<T>&res, T coefa = T(1), T coefb = T(0))
{
  if(coefb != T(0)) {
    res *= coefb;
  } else {
    res = T(0);
  }

  if(c1 == 'N' && c2 == 'N')
    res += coefa*product(a, b);
  else if(c1 == 'T' && c2 == 'N') {
    Matrix<T> aux1(a.rows(),a.cols()); aux1 = transpose(a);
    res += coefa*product(aux1, b);
  } else if(c1 == 'T' && c2 == 'T'){
    Matrix<T> aux1(a.rows(),a.cols()); aux1 = transpose(a);
    Matrix<T> aux2(b.rows(),b.cols()); aux2 = transpose(b);
    res += coefa*product(aux1, aux2);
  } else if(c1 == 'N' && c2 == 'T') {
    Matrix<T> aux2(b.rows(),b.cols()); aux2 = transpose(b);
    res += coefa*product(a, aux2);
  } else if(c1 == 'N' && c2 == 'C') {
    Matrix<T> aux2(b.rows(),b.cols()); aux2 = ctranspose(b);
    res += coefa*product(a, aux2);
  } else if(c1 == 'C' && c2 == 'N') {
    Matrix<T> aux1(a.rows(),a.cols()); aux1 = ctranspose(a);
    res += coefa*product(aux1, b);
  } else if(c1 == 'C' && c2 == 'C') {
    Matrix<T> aux1(a.rows(),a.cols()); aux1 = ctranspose(a);
    Matrix<T> aux2(b.rows(),b.cols()); aux2 = ctranspose(b);
    res += coefa*product(aux1, aux2);
  } else if(c1 == 'T' && c2 == 'C') {
    Matrix<T> aux1(a.rows(),a.cols()); aux1 = transpose(a);
    Matrix<T> aux2(b.rows(),b.cols()); aux2 = ctranspose(b);
    res += coefa*product(aux1, aux2);
  } else if(c1 == 'C' && c2 == 'T') {
    Matrix<T> aux1(a.rows(),a.cols()); aux1 = ctranspose(a);
    Matrix<T> aux2(b.rows(),b.cols()); aux2 = transpose(b);
    res += coefa*product(aux1, aux2);
  }
}

template<class T>
void
matrix_vector_product(const char &c, 
                      const dmtk::Matrix<T>& a, const dmtk::Vector<T> &v, 
		      dmtk::Vector<T>&res, T coefa = T(1), T coefb = T(0))
{
  if(coefb != T(0)) {
    res *= coefb;
  } else {
    res = T(0);
  }

  if(c == 'N'){
    res += coefa*product(a, v);
  }else if(c == 'T'){
    Matrix<T> aux(a.rows(),a.cols()); aux = transpose(a);
    res += coefa*product(aux, v);
  } else if(c == 'C'){
    Matrix<T> aux(a.rows(),a.cols()); aux = ctranspose(a);
    res += coefa*product(aux, v);
  }
}

void
vector_scalar_product(Vector<double> &y, const double& a, Vector<double>&x)
{
  y += a * x;
}

#endif // WITH_LAPACK

} // namespace dmtk

#endif // __DMTK_LA_UTIL_H__
