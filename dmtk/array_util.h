#ifndef __DMTK_ARRAY_UTIL_H__
#define __DMTK_ARRAY_UTIL_H__

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

double
dot_product(const int& _n, 
            const float *v1, const int& stride1, 
            const float *v2, const int& stride2)
{
  DMTK_int n = _n; 
#ifdef WITH_CBLAS
  return cblas_sdot(n, v1, stride1, v2, stride2);  
#else
  return sdot_(n, v1, stride1, v2, stride2);  
#endif
}

double
dot_product(const int& _n, 
            const double *v1, const int& stride1, 
            const double *v2, const int& stride2)
{
  DMTK_int n = _n; 
#ifdef WITH_CBLAS
  return cblas_ddot(n, v1, stride1, v2, stride2);  
#else
  return ddot_(n, v1, stride1, v2, stride2);  
#endif
}

complex<double>
dot_product(const int& _n, 
            const complex<double> *v1, const int& stride1, 
            const complex<double> *v2, const int& stride2)
{
  DMTK_int n = _n; 
  complex<double> sum(0.f,0.f);
  const complex<double> *pi = v1;
  const complex<double> *pj = v2;

//  sum = zdotc_(n, v1, stride1, v2, stride2);  
//  return sum;

  if(stride1 == 1 && stride2 == 1)
    for(int k = 0; k < n; k++)
      sum += std::conj(*pi++)*(*pj++);
  else
    for(int k = 0; k < n; k++){
      sum += std::conj(*pi)*(*pj);
      pi += stride1;
      pj += stride2;
    }

  return sum;
}

void
array_copy(int _n, const double *in, double *out)
{
  DMTK_int n = _n; 
#ifdef WITH_CBLAS
  cblas_dcopy(n, in, 1, out, 1);
#else
  dcopy_(n, in, 1, out, 1);
#endif
}

void
array_copy(int _n, const complex<double> *in, complex<double> *out)
{
  DMTK_int n = _n; 
#ifdef WITH_CBLAS
  cblas_zcopy(n, in, 1, out, 1);
#else
  zcopy_(n, in, 1, out, 1);
#endif
}

template<class T>
void
array_copy(int n, const T* in, T* out)
{
   int m = n%7;
   if(m != 0) {
      for(int i = 0; i < m; i++) out[i] = in[i];
      if(n < 7) return;
   } 
   for(int i = m; i < n; i+= 7){
     out[i] = in[i];
     out[i + 1] = in[i + 1];
     out[i + 2] = in[i + 2];
     out[i + 3] = in[i + 3];
     out[i + 4] = in[i + 4];
     out[i + 5] = in[i + 5];
     out[i + 6] = in[i + 6];
   }
}

#else // !WITH_LAPACK

template<class T>
T
dot_product(const int& n, 
            const T *v1, const int& stride1, const T *v2, const int& stride2)
{
  T sum(0);
  const T *pi = v1;
  const T *pj = v2;

  if(stride1 == 1 && stride2 == 1)
    for(int k = 0; k < n; k++)
      sum += std::conj(*pi++)*(*pj++);
  else
    for(int k = 0; k < n; k++){
      sum += std::conj(*pi)*(*pj);
      pi += stride1;
      pj += stride2;
    }

  return sum;
}

template<class T>
void
array_copy(int n, const T* in, T* out)
{
   int m = n%7;
   if(m != 0) {
      for(int i = 0; i < m; i++) out[i] = in[i];
      if(n < 7) return;
   } 
   for(int i = m; i < n; i+= 7){
     out[i] = in[i];
     out[i + 1] = in[i + 1];
     out[i + 2] = in[i + 2];
     out[i + 3] = in[i + 3];
     out[i + 4] = in[i + 4];
     out[i + 5] = in[i + 5];
     out[i + 6] = in[i + 6];
   }
}

#endif // WITH_LAPACK

template<class T>
void
array_copy(int n, T& in, T& out)
{
   int m = n%7;
   if(m != 0) {
      for(int i = 0; i < m; i++) out[i] = in[i]; 
      if(n < 7) return;
   } 
   for(int i = m; i < n; i+= 7){
     out[i] = in[i];
     out[i + 1] = in[i + 1];
     out[i + 2] = in[i + 2];
     out[i + 3] = in[i + 3];
     out[i + 4] = in[i + 4];
     out[i + 5] = in[i + 5];
     out[i + 6] = in[i + 6];
   }
}

template<class T1, class T2>
void
array_copy2(int n, const T1& in, T2& out)
{
   int m = n%7;
   if(m != 0) {
      for(int i = 0; i < m; i++) out[i] = in[i];
      if(n < 7) return;
   } 
   for(int i = m; i < n; i+= 7){
     out[i] = in[i];
     out[i + 1] = in[i + 1];
     out[i + 2] = in[i + 2];
     out[i + 3] = in[i + 3];
     out[i + 4] = in[i + 4];
     out[i + 5] = in[i + 5];
     out[i + 6] = in[i + 6];
   }
}

template<class T1, class T2>
void
array_copy2(int n, T1& in, T2& out)
{
   int m = n%7;
   if(m != 0) {
      for(int i = 0; i < m; i++) out[i] = in[i];
      if(n < 7) return;
   } 
   for(int i = m; i < n; i+= 7){
     out[i] = in[i];
     out[i + 1] = in[i + 1];
     out[i + 2] = in[i + 2];
     out[i + 3] = in[i + 3];
     out[i + 4] = in[i + 4];
     out[i + 5] = in[i + 5];
     out[i + 6] = in[i + 6];
   }
}

inline double quickran(long & idum)
{
    const int im = 134456;
    const int ia = 8121;
    const int ic = 28411;
    const double scale = 1.0 / im;
    idum = (idum*ia+ic)%im;
    return double(idum) * scale;
}

inline double quickran(long long & idum)
{
    const int im = 134456;
    const int ia = 8121;
    const int ic = 28411;
    const double scale = 1.0 / im;
    idum = (idum*ia+ic)%im;
    return double(idum) * scale;
}

} // namespace dmtk

#endif // __DMTK_ARRAY_UTIL_H__
