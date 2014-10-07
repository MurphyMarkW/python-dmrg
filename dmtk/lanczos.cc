#ifndef __DMTK_LANCZOS_CC__
#define __DMTK_LANCZOS_CC__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "range.h"
#include "vector.h"
#include "matrix.h"
#ifdef WITH_BZIP2
#include "bzip2stream.h"
#endif

using namespace std;

namespace dmtk
{

#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
#define SWAP(a,b) { itemp=(a);(a)=(b);(b)=itemp; }

template<class T>
class AuxIndex
{
  public:
    T v;
    size_t index;

    AuxIndex():v(T()), index(0) {}
    AuxIndex(const T& _v, const size_t &i): v(_v), index(i) {}
    AuxIndex & operator=(const AuxIndex& o) { v = o.v; index = o.index; return *this; } 

    bool operator>(const AuxIndex &o) const
      {
        if(v > o.v) return true;
        if(v == o.v && index > o.index) return true;
        return false;
      }
    bool operator<(const AuxIndex &o) const
      {
        if(v < o.v) return true;
        if(v == o.v && index < o.index) return true;
        return false;
      }
    bool operator>=(const AuxIndex &o) const
      {
        return (v >= o.v && index >= o.index); 
      }
    bool operator<=(const AuxIndex &o) const
      {
        return(v <= o.v && index <= o.index);
      }
    bool operator==(const AuxIndex &o) const
      {
        return (v == o.v);
      }
    bool operator!=(const AuxIndex &o) const
      {
        return (v != o.v);
      }
};


template <class T>
int sort2(size_t n, Vector<T>& arr, bool ascending = true)
{
  size_t j,l;
  T itemp;
  int nswaps = 0;
  bool changed = true;
  T a;

  if(n == 0){
    cout << "*** WARNING: sort width (n == 0) \n";
    return -1;
  }

  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  while(changed){
    changed = false;
    for(l = 0; l < n-1; l++){
      if(aux[l] > aux[l+1]){
        SWAP(aux[l],aux[l+1]); 
        nswaps++;
        changed = true;
      }
    }
  }

  if(!ascending && n > 1)
    for(j = 0; j <= n/2-1; j++)
      if(aux[j] != aux[n-j-1]) {
        SWAP(aux[j],aux[n-j-1]);
        nswaps++;
      }

  for(int i = 0; i < n; i++) arr[i] = aux[i].v;  
  return nswaps;
}

template <class T, class A>
int indexx2(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
  size_t j,l,itemp;
  bool changed = true;
  int nswaps = 0;

  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  for (j=0;j<n;j++) indx[j]=j;
  while(changed){
    changed = false;
    for(l = 0; l < n-1; l++){
      if(aux[indx[l]] > aux[indx[l+1]]){
        SWAP(indx[l],indx[l+1]); 
        nswaps++;
        changed = true;
      }
    }
  }

  if(!ascending && n > 1)
    for(j = 0; j <= n/2-1; j++){
      SWAP(indx[j],indx[n-j-1]);
      nswaps++;
    }

  return nswaps;
}

// FROM NR
/*
template <class T>
void sort(size_t n, Vector<T>& arr, bool ascending = true)
{
  if (n < 2) return;

  int l, j, ir, i;
  T itemp;
  AuxIndex<T> rra;
                                                                              
  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  ra = &aux[0]-1; 

  l = n / 2 + 1;
  ir = n;
  for (;;)
      {
      if (l > 1)
          rra = ra[--l];
      else
          {
          rra = ra[ir];
          ra[ir] = ra[1];
          if (--ir == 1)
              {
              ra[1] = rra;
              break;
              }
          }
      i = l;
      j = l + l;
      while (j <= ir)
          {
          if (j < ir && ra[j] < ra[j + 1])
              ++j;
          if (rra < ra[j])
              {
              ra[i] = ra[j];
              j += (i = j);
              }
          else
              j = ir + 1;
          }
      ra[i] = rra;
      }

  ra = &aux[0]; 
  for(i = 0; i < n; i++) { arr[i] = ra->v; ra++; } 

  if(!ascending && n > 1)
    for(j = 0; j <= n/2-1; j++)
      if(arr[j] != arr[n-j-1]) SWAP(arr[j],arr[n-j-1]);
}

template<class T, class A>
void indexx(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
    if (n == 1) { indx[0] = 0; return; }

    int l, j, ir, i;
    AuxIndex<T> rra;
    size_t rrb;
    size_t itemp;

    Vector<AuxIndex<T> > aux(n);
    AuxIndex<T> *ra = &aux[0];
    for(i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

    ra = &aux[0]-1; 
    size_t *rb = indx.array()-1;

    for(i = 1; i <= n; i++) rb[i] = i-1;         
                                                                       
    l = (n >> 1) + 1;
    ir = n;
    for (;;)
        {
        if (l > 1)
            {
            rra = ra[--l];
            rrb = rb[l];
            }
        else
            {
            rra = ra[ir];
            rrb = rb[ir];
            ra[ir] = ra[1];
            rb[ir] = rb[1];
            if (--ir == 1)
                {
                ra[1] = rra;
                rb[1] = rrb;
                break;
                }
            }
        i = l;
        j = l << 1;
        while (j <= ir)
            {
            if (j < ir && ra[j] < ra[j + 1])
                ++j;
            if (rra < ra[j])
                {
                ra[i] = ra[j];
                rb[i] = rb[j];
                j += (i = j);
                }
            else
                j = ir + 1;
            }
        ra[i] = rra;
        rb[i] = rrb;
        }

  if(!ascending && n > 1)
    for(j = 0; j <= n/2-1; j++)
      SWAP(indx[j],indx[n-j-1]);
}
*/
#undef M
#undef NSTACK
// END NR

template <class T>
void sort(size_t n, Vector<T>& arr, bool ascending = true)
{
  if (n < 2) return;

  T itemp;

  std::vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; }

  std::sort(aux.begin(), aux.end());

  ra = &aux[0];
  for(int i = 0; i < n; i++, ra++) arr[i] = ra->v;

  if(!ascending && n > 1)
    for(int j = 0; j <= n/2-1; j++)
      if(arr[j] != arr[n-j-1]) SWAP(arr[j],arr[n-j-1]);
}

template<class T, class A>
void indexx(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true){
  if (n == 1) { indx[0] = 0; return; }

  size_t itemp;

  std::vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; }

  std::sort(aux.begin(), aux.end());

  ra = &aux[0];
  for(int i = 0; i < n; i++, ra++) indx[i] = ra->index;

  if(!ascending && n > 1)
    for(int j = 0; j <= n/2-1; j++)
      SWAP(indx[j],indx[n-j-1]);
}

#ifdef WITH_LAPACK
#ifndef WITH_NR
void tqli(dmtk::Vector<double>& d, dmtk::Vector<double>& e,
          int _n, dmtk::Matrix<double>& z, bool calc_vectors = false)
{
  char calc = calc_vectors ? 'V' : 'N';
  DMTK_int n = _n;
  DMTK_int lwork = calc_vectors ? 1 + 4*n + n*n : 1;
  Vector<double> work(lwork);
  DMTK_int liwork = calc_vectors ? 3 + 5*n : 1;
  Vector<DMTK_int> iwork(liwork);
  Vector<double> auxd(n);
  Vector<double> auxe(n);
  Matrix<double> auxz(n,n);
  auxz(Range(0,n-1),Range(0,n-1)) = z(Range(1,n),Range(1,n));
  auxd(Range(0,n-1)) = d(Range(1,n));
  auxe(Range(0,n-1)) = e(Range(2,n+1));

  DMTK_int info = 0;
  dstevd_(calc,n,auxd.array(),auxe.array(),auxz.array(),n,work.array(),lwork,iwork.array(),liwork,info);
  if(info != 0) cout << "ERROR in dstevd " << info << endl;

  z(Range(1,n),Range(1,n)) = auxz(Range(0,n-1),Range(0,n-1));
  d(Range(1,n)) = auxd(Range(0,n-1));
  e(Range(1,n)) = auxe(Range(0,n-1));
}

#else // !WITH_LAPACK

/* (C) Copr. 1986-92 Numerical Recipes Software . */
void tred2(dmtk::Matrix<double> &a, int n,
           dmtk::Vector<double> &d, dmtk::Vector<double> &e,
           bool calc_vectors = false)
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0f;
    if (l > 1) {
      for (k=1;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0f)
        e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0f ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0f;
  e[1]=0.0f;
  /* Contents of this loop can be omitted if eigenvectors not
      wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l=i-1;
    if(calc_vectors){
      if (d[i]) {
        for (j=1;j<=l;j++) {
          g=0.0f;
          for (k=1;k<=l;k++)
            g += a[i][k]*a[k][j];
          for (k=1;k<=l;k++)
            a[k][j] -= g*a[k][i];
        }
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0f;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}

/* (C) Copr. 1986-92 Numerical Recipes Software . */
void tqli(dmtk::Vector<double>& d, dmtk::Vector<double>& e, 
          int n, dmtk::Matrix<double>& z, bool calc_vectors = false)
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0f;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
        if (iter++ == 30)  cout << "Too many iterations in TQLI" << endl;
        g=(d[l+1]-d[l])/(2.0f*e[l]);
        r=sqrt((g*g)+1.0f);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0f;
        p=0.0f;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          if (fabs(f) >= fabs(g)) {
            c=g/f;
            r=sqrt((c*c)+1.0f);
            e[i+1]=f*r;
            c *= (s=1.0f/r);
          } else {
            s=f/g;
            r=sqrt((s*s)+1.0f);
            e[i+1]=g*r;
            s *= (c=1.0f/r);
          }
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0f*c*b;
          p=s*r;
          d[i+1]=g+p;
          g=c*r-b;
          /* Next loop can be omitted if eigenvectors not wanted */
          if(calc_vectors){
            for (k=1;k<=n;k++) {
              f=z[k][i+1];
              z[k][i+1]=s*z[k][i]+c*f;
              z[k][i]=c*z[k][i]-s*f;
            }
          }
        }
        d[l]=d[l]-p;
        e[l]=g;
        e[m]=0.0f;
      }
    } while (m != l);
  }
}

#endif // WITH_NR
#endif //WITH_LAPACK

//////////////////////////////////////////////////////////////////////

template<class T, class A, class B>
void lanczos(A& m,
             Vector<B>& gs, const B& seed, 
             Vector<double>& ener, int nvectors,
             dmtk::Vector<double>& a, dmtk::Vector<double>& b, 
             int &maxiter, double tol,
             bool use_seed, bool calc_gs, const char* vectors_file, bool force_maxiter = false)
{  
  B x1(gs(0)), x2(gs(0)), aux(gs(0));
  dmtk::Vector<double> d(1000), e(1000), xc(1000);
  dmtk::Matrix<double> z(1000,1000);
  uint col, iter = 0;
  double eini, e0;
  int control_max = maxiter;

  if(maxiter == -1) force_maxiter = false;

  if(control_max == 0) { gs[0] = T(1); maxiter = 1; return; }
#ifndef USE_LANCZOS_VECTORS
  ofstream real_file(vectors_file,std::ios::out|std::ios::binary);
  if(!real_file) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
#ifdef WITH_BZIP2
  bzip2_stream::bzip2_ostream outputfile(real_file);
#else
  ofstream &outputfile = real_file;
#endif
#endif
 
  e0 = 10000.0f;
  x1 = x2 = aux = T(0);
  maxiter = 0;

  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    if(use_seed){
      x1 = seed;
    } else { 
#ifdef USE_LANCZOS_BASIC
      x1 = 1.;
//      x1 = 0;
//      x1[0] = 1.;
#else
      x1.randomize();
#endif
    } 

    iter = 0;

    b(0) = sqrt(real(dot_product(x1,x1)));
    x1 = x1 / T(b(0));
    x2 = T(0);
    b(0) = 1.;


    uint nmax = std::min(999, (int)gs(0).size());
    for(iter = 1; iter <= nmax; iter++){ // Lanczos iteration
      eini = e0;

      if(b(iter - 1) != 0.){
        aux = x1;
        x1 = -T(b(iter-1)) * x2;
        x2 = aux / T(b(iter-1));
      }

//---------------------------------------------------------------
//******* REORTHOGONALIZATION
/*
      outputfile.close();
      ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
      if(!inputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";

      B aux1 = x2, aux2 = x2;
      T overlap;
      if(iter >= 1){ 
         for(int i = 1; i < iter; i++){
           aux2.read(inputfile);
           overlap = product(aux1,aux2);
           x2 -= overlap*aux2;
         }
         x2 /= product(x2,x2);

//         inputfile.close();
//         ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
//         for(int i = 1; i < iter; i++){
//           aux2.read(inputfile);
//           overlap = product(x2,aux2);
//           cout << i << " " << overlap << endl;
//         }

      }

      inputfile.close();
      ofstream outputfile(vectors_file,std::ios::out|std::ios::binary|std::ios::app);
      if(!outputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
*/
//---------------------------------------------------------------
      aux = product(m,x2);
      x1 = x1 + aux;
      T xa = dot_product(x1, x2);
      a(iter) = real(xa);
      x1 = x1 - x2*T(a(iter));
      b(iter) = sqrt(real(dot_product(x1, x1)));
   
//      cout << "Iter=" << iter << " " << a(iter) << " " << b(iter) << endl;

#ifndef USE_LANCZOS_VECTORS
      if(calc_gs) x2.write(outputfile);
#endif

      if(maxiter > 0) {  // We are building the ground state;

#ifdef USE_LANCZOS_VECTORS
        gs = gs + xc(iter) * x2;
        if(iter == maxiter) return; 
#endif
      } else{  // we are calculating energy

        if(iter >= 2){
          d(Range(1,iter)) = a(Range(1,iter));
          e(Range(2,iter+1)) = b(Range(1,iter));

          // call tqli without eigenvectors
          tqli(d, e, iter, z, false);
  
#ifdef WITH_LAPACK
          e0 = d(1);
          col = 1;
#endif
#ifdef WITH_NR
          e0 = 10000.0f;
          for(int j = 1; j <= iter; j++){
            if(d(j) < e0) {
              e0 = d(j);
              col = j;
            }
          }
#endif

          cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << endl;
          if((force_maxiter && iter >= control_max) || (iter >= gs(0).size()-1 || iter == 999 || fabs(b(iter)) < tol) || (force_maxiter == false && fabs(eini-e0) <= tol))
          { 
             // converged
             if(fabs(b(iter)) < tol) cout << "b(iter) < tol " << b(iter) << endl;
             cout << setprecision(12) << "E0 = " << e0 << endl;
             maxiter = iter;
             if(!calc_gs) return; // We return with ground states energy
             z = I<double>(); //identity;
             d(Range(1,iter)) = a(Range(1,iter));
             e(Range(2,iter+1)) = b(Range(1,iter));
             // call tqli with eigenvectors
             tqli(d, e, iter, z, true);
             xc = z(col, Range(0,iter));

#ifdef USE_LANCZOS_VECTORS
             break; // Exit Lanczos iteration. Re-start for ground state.
#endif
             outputfile.flush();
             outputfile.close();
         
             Vector<size_t> indx(maxiter); 
             d(Range(0,maxiter-1)) = d(Range(1,maxiter));
             indexx<double,Vector<double> >(maxiter, d, indx, true);
             for(int n = 0; n < maxiter; n++)
               cout << "LEVEL " << n << " " << " " << d(indx(n)) << endl;

             for(int n = 0; n < nvectors; n++){
               ifstream real_file(vectors_file,std::ios::in|std::ios::binary);
               if(!real_file) cerr << "*** ERROR 2: Lanczos could not open vectors.dat\n";
#ifdef WITH_BZIP2
               bzip2_stream::bzip2_istream inputfile(real_file);
#else
               ifstream &inputfile = real_file;
#endif
               gs(n) = T(0);
               ener(n) = d(indx(n));
//               cout << n << " " << ener(n) << endl;
               for(int i = 1; i <= maxiter; i++){
                 x2.read(inputfile);
                 gs(n) += T(z(indx(n)+1,i)) * x2;
               }
               inputfile.close();

             }

             return;
          }
        } // diagonalization of tridiagonal matrix

      } 
    } // Lanczos iteration
  } // main iteration

}

//////////////////////////////////////////////////////////////////////

template<class T, class A, class B>
void lanczos(A& m,
             B& gs, const B& seed, 
             double& ener, 
             dmtk::Vector<double>& a, dmtk::Vector<double>& b, 
             int &maxiter, double tol,
             bool use_seed, bool calc_gs, const char* vectors_file, bool force_maxiter = false)
{  
  B x1(gs), x2(gs), aux(gs);
  dmtk::Vector<double> d(1000), e(1000), xc(1000);
  dmtk::Matrix<double> z(1000,1000);
  uint col, iter = 0;
  double eini, e0;
  int control_max = maxiter;

  if(maxiter == -1) force_maxiter = false;

  if(control_max == 0) { gs = T(1); maxiter = 1; return; }
#ifndef USE_LANCZOS_VECTORS
  ofstream real_file(vectors_file,std::ios::out|std::ios::binary);
  if(!real_file) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
#ifdef WITH_BZIP2
  bzip2_stream::bzip2_ostream outputfile(real_file);
#else
  ofstream &outputfile = real_file;
#endif
#endif
  e0 = 10000.0f;
  x1 = x2 = aux = T(0);

  maxiter = 0;
  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    if(use_seed){
      x1 = seed;
    } else { 
#ifdef USE_LANCZOS_BASIC
      x1 = 1.;
//      x1 = 0;
//      x1[0] = 1.;
#else
      x1.randomize();
#endif
    } 

    iter = 0;

    b(0) = sqrt(real(dot_product(x1,x1)));
    x1 = x1 / T(b(0));
    x2 = T(0);
    b(0) = 1.;

    uint nmax = std::min(999, (int)gs.size());
    for(iter = 1; iter <= nmax; iter++){ // Lanczos iteration
      eini = e0;

      if(b(iter - 1) != 0.){
        aux = x1;
        x1 = -T(b(iter-1)) * x2;
        x2 = aux / T(b(iter-1));
      }

      aux = T(0);
      aux = product(m,x2);
      x1 = x1 + aux;
      a(iter) = real(dot_product(x1, x2));
      x1 = x1 - x2*T(a(iter));
      b(iter) = sqrt(real(dot_product(x1, x1)));
   
//      cout << "Iter=" << iter << " " << a(iter) << " " << b(iter) << endl;

#ifndef USE_LANCZOS_VECTORS
      if(calc_gs) x2.write(outputfile);
#endif

      if(maxiter > 0) {  // We are building the ground state;

#ifdef USE_LANCZOS_VECTORS
        gs = gs + xc(iter) * x2;
        if(iter == maxiter) return; 
#endif

      } else{  // we are calculating energy

        if(iter >= 2){

          d(Range(1,iter)) = a(Range(1,iter));
          e(Range(2,iter+1)) = b(Range(1,iter));
          // call tqli without eigenvectors
          tqli(d, e, iter, z, false);

#ifdef WITH_LAPACK
          e0 = d(1);
          col = 1;
#endif
#ifdef WITH_NR
          e0 = 10000.0f;
          for(int j = 1; j <= iter; j++){
            if(d(j) < e0) {
              e0 = d(j);
              col = j;
            }
          }
#endif

          cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << endl;
          if((force_maxiter && iter >= control_max) || (iter >= gs.size()-1 || iter == 999 || fabs(b(iter)) < tol) || (!force_maxiter  && fabs(eini-e0) <= tol))
          { 
             // converged
             cout << setprecision(12) << "E0 = " << e0 << endl;
             maxiter = iter;
             if(!calc_gs) return; // We return with ground states energy
             z = I<double>(); //identity;
             d(Range(1,iter)) = a(Range(1,iter));
             e(Range(2,iter+1)) = b(Range(1,iter));
             // call tqli with eigenvectors
             tqli(d, e, iter, z, true);
             xc = z(col, Range(0,iter));

#ifdef USE_LANCZOS_VECTORS
             break; // Exit Lanczos iteration. Re-start for ground state.
#endif

             outputfile.flush();
             outputfile.close();

             Vector<size_t> indx(maxiter); 
             d(Range(0,maxiter-1)) = d(Range(1,maxiter));
             indexx<double,Vector<double> >(maxiter, d, indx, true);
             for(int n = 0; n < maxiter; n++)
               cout << "LEVEL " << n << " " << " " << d(indx(n)) << endl;
         
             ifstream real_file(vectors_file,std::ios::in|std::ios::binary);
             if(!real_file) cerr << "*** ERROR 2: Lanczos could not open vectors.dat\n";
#ifdef WITH_BZIP2
             bzip2_stream::bzip2_istream inputfile(real_file);
#else
             ifstream &inputfile = real_file;
#endif
             gs = T(0);
             ener = e0;
             for(int i = 1; i <= maxiter; i++){
               x2.read(inputfile);
               gs += T(xc(i)) * x2;
             }

             inputfile.close();

             return;
          }
        } // diagonalization of tridiagonal matrix

      } 
    } // Lanczos iteration
  } // main iteration

}

template<class T, class A, class B>
void lanczos_nonsymmetric(A& m,
             B& gs, const B& seed, 
             double& ener, 
             dmtk::Vector<double>& a, dmtk::Vector<double>& b, dmtk::Vector<double> g,
             int &maxiter, double tol,
             bool use_seed, bool calc_gs, const char* vectors_file, bool force_maxiter = false)
{  
  B r(gs), v(gs), s(gs), w(gs), aux(gs);
  dmtk::Matrix<double> h(1000,1000);
  Vector<double> xc(1000), d(1000);
  uint col, iter = 0;
  double eini, e0;
  int control_max = maxiter;

  if(maxiter == -1) force_maxiter = false;

  if(control_max == 0) { gs = T(1); maxiter = 1; return; }
#ifndef USE_LANCZOS_VECTORS
  ofstream real_file(vectors_file,std::ios::out|std::ios::binary);
  if(!real_file) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
#ifdef WITH_BZIP2
  bzip2_stream::bzip2_ostream outputfile(real_file);
#else
  ofstream &outputfile = real_file;
#endif
#endif
  e0 = 10000.0f;
  r = T(0);
  s = T(0);
  maxiter = 0;
  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    g = 0.0f;
    if(use_seed){
      v = seed;
      r = v;
      w = seed;
      s = w;
    } else { 
#ifdef USE_LANCZOS_BASIC
      x1 = 1.;
//      x1 = 0;
//      x1[0] = 1.;
#else
      r.randomize();
      s.randomize();
#endif
    } 

    iter = 0;

    uint nmax = std::min(999, (int)gs.size());
    double rho = 1.;
    double eta = 1.;
    double old_delta = 1.;
    h = T(0);

    b(0) = sqrt(product(r,r));
    g(0) = product(s,r)/b(iter);
    v = 0.;
    w = 0.;
    for(iter = 1; iter <= nmax; iter++){ // Lanczos iteration
      eini = e0;

      aux = v;
      v = r/b(iter-1);
      r = - g(iter-1)*aux; 

      aux = w;
      w = s/g(iter-1);
      s = -b(iter-1)*aux;
      
      aux = product(m,v);
      a(iter) = product(w,aux);

      r += aux - a(iter)*v;

      aux = product_hc(m,w);
      s += aux - a(iter)*w;
// cout << "HOLA OVERLAP " << product(v,w) << " " << a(iter) << " " << product(v,aux) << endl;
      b(iter) = sqrt(product(r,r));
      g(iter) = product(s,r)/b(iter);

      h(iter-1,iter-1) = a(iter);
      h(iter, iter-1) = b(iter);
      h(iter-1, iter) = g(iter);
      
//      cout << "Iter=" << iter << " " << a(iter) << " " << b(iter) << " " << g(iter) << endl;

#ifndef USE_LANCZOS_VECTORS
      if(calc_gs) {
        aux = v; 
        aux.write(outputfile);
//        w.write(outputfile);
      }
#endif

      if(maxiter > 0) {  // We are building the ground state;

#ifdef USE_LANCZOS_VECTORS
        gs = gs + xc(iter-1) * v; //r/b(iter);
        if(iter == maxiter) return; 
#endif

      } else{  // we are calculating energy

        if(iter >= 2){
          Vector<double> wr(iter);
          Vector<double> wi(iter);
          Matrix<double> realh(iter,iter);
          for(int i = 0; i < iter-1; i++){
            realh(i,i) = h(i,i);
            realh(i+1,i) = h(i+1,i); 
            realh(i,i+1) = h(i,i+1); 
          }
          realh(iter-1,iter-1) = h(iter-1,iter-1);

          DMTK_int info;
          DMTK_int ilo, ihi;
          Vector<double> scale(iter);
          Vector<double> work(2*iter);
          Vector<DMTK_int> iwork((2*iter-2 >1) ? 2*iter-2 : 1);
          Vector<double> rconde(iter), rcondv(iter);
          dmtk::Matrix<double> zl(iter,iter),zr(iter,iter);
          double abnrm;
          dgeevx_('N','N','N','N',iter, realh.array(), iter, wr.array(), wi.array(), zl.array(), iter, zr.array(), iter, ilo, ihi, scale.array(), abnrm, rconde.array(), rcondv.array(), work.array(), work.size(), iwork.array(), info);

          Vector<size_t> indx(iter); 
          d(Range(0,iter-1)) = wr(Range(0,iter-1));
          indexx<double,Vector<double> >(iter, d, indx, true);
          e0 = d(indx(0));
          col = indx(0);
         
          cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << " " << wr(indx(0)) << " " << wi(indx(0)) << endl;

          if((force_maxiter && iter >= control_max) || (iter >= gs.size()-1 || iter == 999 || fabs(b(iter)) < tol) || (!force_maxiter  && fabs(eini-e0) <= tol))
          { 
             // converged
             cout << setprecision(12) << "E0 = " << e0 << endl;
             maxiter = iter;
             if(!calc_gs) return; // We return with ground states energy

             realh = 0.; 
             for(int i = 0; i < iter-1; i++){
               realh(i,i) = h(i,i);
               realh(i+1,i) = h(i+1,i); 
               realh(i,i+1) = h(i,i+1); 
             }
             realh(iter-1,iter-1) = h(iter-1,iter-1);
             work.resize(iter*(iter+6));
             dgeevx_('N','N','V','N',iter, realh.array(), iter, wr.array(), wi.array(), zl.array(), iter, zr.array(), iter, ilo, ihi, scale.array(), abnrm, rconde.array(), rcondv.array(), work.array(), work.size(), iwork.array(), info);

             d(Range(0,iter-1)) = wr(Range(0,iter-1));
             indexx<double,Vector<double> >(iter, d, indx, true);
             col = indx(0);

             xc = zr(col, Range(0,iter-1));

#ifdef USE_LANCZOS_VECTORS
             break; // Exit Lanczos iteration. Re-start for ground state.
#endif
             outputfile.flush();
             outputfile.close();

             for(int n = 0; n < maxiter; n++)
               cout << "LEVEL " << n << " " << " " << d(indx(n)) << endl;
         
             ifstream real_file(vectors_file,std::ios::in|std::ios::binary);
             if(!real_file) cerr << "*** ERROR 2: Lanczos could not open vectors.dat\n";

#ifdef WITH_BZIP2
             bzip2_stream::bzip2_istream inputfile(real_file);
#else
             ifstream &inputfile = real_file;
#endif
             gs = T(0);
             ener = e0;
             aux = 0.;
             for(int i = 0; i < maxiter; i++){
               r.read(inputfile);
               gs += T(xc(i)) * r;
             }
             gs /= sqrt(product(gs,gs));
             aux = product(m,gs);
             cout << "HOLA ENERGY " << product(gs,gs) << " " << product(gs,aux) << endl;

             inputfile.close();
             return;
          }
        } // diagonalization of tridiagonal matrix

      } 
    } // Lanczos iteration
  } // main iteration

}

template<class T, class A, class B>
void power_method(A& m,
                  B& gs, const B& seed, 
                  double& ener, 
                  int &maxiter, double tol,
                  bool use_seed, bool force_maxiter = false)
{  
  B v(gs), w(gs);
  int iter = 0;
  double eps = 1.;

  if(use_seed)
    v = seed;
  else
    v.randomize();

  v /= sqrt(product(v,v));

  while(true){
    iter++;

    w = product(m,v);
    double e = real(product(v,w));
    w *= T(eps);
    w += v*T(1.-eps);
//    w *= T(-1.);
    v = w*T(1./sqrt(product(w,w)));

    cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e << endl;

    if((force_maxiter && iter == maxiter) || fabs(e-ener) < tol) {
      e = ener;
      gs = v;
      cout << setprecision(12) << "E0 = " << e << endl;
      break;
    } 
    ener = e;
  }
  return;
}

template<class T, class A, class B>
void power_method_modified(A& m,
                  B& gs, const B& seed, 
                  double& ener, 
                  int &maxiter, double tol,
                  bool use_seed, bool force_maxiter = false)
{  
  B v1(gs), v2(gs), w(gs), aux(gs);
  int iter = 0;
  double a, b, c, e = 10000;

  if(use_seed)
    v1 = seed;
  else
    v1.randomize();
  
  v1 /= sqrt(product(v1,v1));

  while(true){
    iter++;

    aux = v1;
    w = product(m,v1);
    a = product(w,w);
    v2 = w - product(w,v1)*v1;
    v2 /= sqrt(product(v2,v2));

    v1 = product(m,v2);
    c = product(v1,w);
    b = product(v1,v1);

    e = 0.5*((a+b)-sqrt((a-b)*(a-b)+4*c*c));
    double x = (e-a)/c;
    v1 = (aux+v2*x)*(1./sqrt(1.+x*x)); 
    
    cout << setprecision(12) << "Iter = " << iter << " " << "  Ener = " << e << endl;

    if((force_maxiter && iter == maxiter) || fabs(e-ener) < tol) {
      gs = v1;
      w = product(m,gs);
      e = product(gs,w);
      cout << setprecision(12) << "E0 = " << e << endl;
      break;
    } 
    ener = e;
  }
  return;
}

template<class T, class A, class B>
void lanczos_modified(A& m,
                  B& gs, const B& seed, 
                  double& ener, 
                  int &maxiter, double tol,
                  bool use_seed, bool force_maxiter = false)
{  
  B v1(gs), v2(gs), w(gs), aux(gs);
  int iter = 0;
  double a, b, c, e = 10000;

  if(use_seed)
    v1 = seed;
  else
    v1.randomize();
  
  v1 /= sqrt(product(v1,v1));

  while(true){
    iter++;

    aux = v1;
    w = product(m,v1);
    a = product(v1,w);
    v2 = w - a*v1;
    v2 /= sqrt(product(v2,v2));

    v1 = product(m,v2);
    b = product(v1,v2);
    double c1 = product(aux,v1);
    double c2 = product(v2,w);

    e = 0.5*((a+b)-sqrt((a-b)*(a-b)+4*c1*c2));
    double x = (e-a)/c1;
    v1 = (aux+v2*x)*(1./sqrt(1.+x*x)); 
//    cout << "HAMI " << a << " " << b << " " <<  c1 << " " << c2 << endl;
    cout << setprecision(12) << "Iter = " << iter << " " << "  Ener = " << e << endl;

    if((force_maxiter && iter == maxiter) || fabs(e-ener) < tol) {
      gs = v1;
      gs /= sqrt(product(gs,gs));
      w = product(m,gs);
      e = product(gs,w);
      cout << setprecision(12) << "E0 = " << e << endl;
      break;
    } 
    ener = e;
  }
  return;
}

template<class T, class A, class B>
void lanczos_unsymmetric_modified(A& m,
                  B& gs, B& gs_left, const B& seed, 
                  double& ener, double &eneri, 
                  int &maxiter, double tol,
                  bool use_seed, bool force_maxiter = false, int mask_hc = (MASK_PRODUCT_DEFAULT|MASK_PRODUCT_HC))
{  
  B v1(gs), v2(gs), w(gs), aux(gs);
  int iter = 0;
  T a, b, c;
  double e = 10000, ei = 10000;

  if(use_seed)
    v1 = seed;
  else
    v1.randomize();
  
  v1 /= sqrt(product(v1,v1));

  while(true){
    iter++;

    aux = v1;
    w = product(m,v1);
    a = product(v1,w);
    v2 = w - a*v1;
    v2 /= sqrt(product(v2,v2));

    v1 = product(m,v2);
    b = product(v1,v2);
    T c1 = product(aux,v1);
    T c2 = product(v2,w);

    ei = 0.;
    T arg = (a-b)*(a-b)+T(4.)*c1*c2;
    T xe = T(0.5)*((a+b)-sqrt(arg)); 
    e = real(xe);
    ei = abs(imag(xe));
    xe = e + ei;
    
    T x = (xe-a)/c1;
    v1 = (aux+v2*x);
    v1 /= sqrt(product(v1,v1)); 
//cout << "HOLA " << product(aux,aux) << " " << product(v2,v2) << " " << product(v1,v1) << endl;
//cout << "HOLA HAMI " << a << " " << b << " " <<  c1 << " " << c2 << endl;
    cout << setprecision(12) << "Modified Iter = " << iter << " " << "  Ener = " << e << " " << ei << endl;

    if((force_maxiter && iter == maxiter) || fabs(e-ener+ei-eneri) < tol) {
      gs = v1;
      x = (xe-a)/c2;
      v1 = (aux+v2*x);
      v1 /= sqrt(product(v1,v1)); 
      gs_left = v1; 
      cout << setprecision(12) << "Modified E0 = " << e << " " << ei << endl;
      break;
    } 
    ener = e;
    eneri = ei;
  }
  return;
}

template<class T, class A, class B>
Matrix<T>
verify_hamiltonian(A& m, const B& gs)
{
  Matrix<T> h(gs.size(),gs.size());
  B aux(gs);

  cout << "***** VERIFYING HAMILTONIAN *****\n";
  for(int col = 0; col < gs.size(); col++){
    aux = T(0);
    aux(col) = T(1);
    h.column(col) = product(m,aux);
//    for(int row = 0; row < gs.size(); row++){
//      cout << col << " " << row << " " << h.column(col)(row) << endl;
//    }
  }
  return h;
}


#undef SWAP
#undef SIGN

} // namespace dmtk

#endif // __DMTK_LANCZOS_CC__
