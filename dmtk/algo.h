#ifndef __DMTK_ALGO_H__
#define __DMTK_ALGO_H__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"

using namespace std;

namespace dmtk
{

#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
#define SWAP(a,b) { itemp=(a);(a)=(b);(b)=itemp; }


template <class T>
void sort(size_t n, Vector<T>& arr, bool ascending = true)
{
  size_t j,l;
  T itemp;
  bool changed = true;
  T a;

  if(n == 0){
    cout << "*** WARNING: sort width (n == 0) \n";
    return;
  }

  while(changed){
    changed = false;
    for(l = 0; l < n-1; l++){
      if(arr[l] > arr[l+1]){
        SWAP(arr[l],arr[l+1]); 
        changed = true;
      }
    }
  }

  if(!ascending)
    for(j = 0; j <= n/2-1; j++)
      if(arr[j] != arr[n-j-1]) SWAP(arr[j],arr[n-j-1]);
}


// FROM NR
#define M 7
#define NSTACK 50
template <class T, class A>
void indexx2(size_t n, const A& arr, Vector<size_t>& indx, bool ascending);

template<class T>
class AuxIndex
{
  public:
    T v;
    size_t index;

    AuxIndex():v(T()), index(0) {}
    AuxIndex(const T& _v, const size_t &i): v(_v), index(i) {}
    AuxIndex & operator=(const AuxIndex& o) { v = o.v; index = o.index; return *this; } 

    bool operator>(const AuxIndex &o)
      {
        if(v > o.v) return true;
        if(v == o.v && index > o.index) return true;
        return false;
      }
    bool operator<(const AuxIndex &o)
      {
        if(v < o.v) return true;
        if(v == o.v && index < o.index) return true;
        return false;
      }
    bool operator>=(const AuxIndex &o)
      {
        return (v >= o.v && index >= o.index); 
      }
    bool operator<=(const AuxIndex &o)
      {
        return(v <= o.v && index <= o.index);
      }
    bool operator==(const AuxIndex &o)
      {
        return (v == o.v);
      }
    bool operator!=(const AuxIndex &o)
      {
        return (v != o.v);
      }
};

template<class T, class A>
void indexx(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
/*
    indexx2<T,A>(n, arr, indx, ascending);
    return; 
*/

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
                                                                       
    if (n < 2)
        return;
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

  if(!ascending)
    for(j = 0; j <= n/2-1; j++)
      SWAP(indx[j],indx[n-j-1]);
}

#undef M
#undef NSTACK
// END NR

template <class T, class A>
void indexx2(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
  size_t j,l,itemp;
  bool changed = true;
  T a;

  for (j=0;j<n;j++) indx[j]=j;
  while(changed){
    changed = false;
    for(l = 0; l < n-1; l++){
      if(arr[indx[l]] > arr[indx[l+1]]){
        SWAP(indx[l],indx[l+1]); 
        changed = true;
      }
    }
  }

  if(!ascending)
    for(j = 0; j <= n/2-1; j++)
      SWAP(indx[j],indx[n-j-1]);

}

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
/* (C) Copr. 1986-92 Numerical Recipes Software . */

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
  dmtk::Vector<double> d(101), e(101), xc(101);
  dmtk::Matrix<double> z(101,101);
  uint i, j, col, iter = 0;
  double eini, e0;
  int control_max = maxiter;
  maxiter = 0;

  ofstream outputfile(vectors_file,std::ios::out|std::ios::binary);
  if(!outputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
       
 
  e0 = 10000.0f;
  x1 = x2 = aux = T(0);

  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    if(use_seed){
      x1 = seed;
    } else { 
      srand(9872121); // arbitrary seed
      for(i = 0; i < x1.size(); i++){
         x1(i) = random() * 1.e-10 - random() * 1.e-10; 
//         x1(i) = 1.;
      }
    } 

    iter = 0;

    b(0) = sqrt(real(product(x1,x1)));
    x1 = x1 / T(b(0));
    x2 = T(0);
    b(0) = 1.;


    uint nmax = std::min(99, (int)gs(0).size()-1);
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
      T xa = product(x1, x2);
      a(iter) = real(xa);
      x1 = x1 - x2*T(a(iter));
      b(iter) = sqrt(real(product(x1, x1)));
   
//      cout << "Iter=" << iter << " " << xa << " " << a(iter) << " " << b(iter) << endl;

      if(calc_gs) x2.write(outputfile);

      if(maxiter > 0) {  // We are building the ground state;

//        gs = gs + xc(iter) * x2;
//        if(iter == maxiter) return; 

      } else{  // we are calculating energy

        if(iter >= 2){

          d(Range(1,iter)) = a(Range(1,iter));
          e(Range(2,iter+1)) = b(Range(1,iter));
          // call tqli without eigenvectors
          tqli(d, e, iter, z, false);

  
          e0 = 10000.0f;
          for(j = 1; j <= iter; j++){
            if(d(j) < e0) {
              e0 = d(j);
              col = j;
            }
          }

          cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << endl;
          if((force_maxiter && iter == control_max) || (!force_maxiter && (iter == gs(0).size()-1 || iter == 99 || fabs(b(iter)) < tol || fabs(eini-e0) <= tol || iter == control_max)))
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

             outputfile.flush();
             outputfile.close();

//             break; // Exit Lanczos iteration. Re-start for ground state.

         
             Vector<size_t> indx(maxiter); 
             d(Range(0,maxiter-1)) = d(Range(1,maxiter));
             indexx<double,Vector<double> >(maxiter, d, indx, true);
//             for(int n = 0; n < maxiter; n++)
//               cout << n << " " << indx(n) << " " << d(n) << " " << d(indx(n)) << endl;

             for(int n = 0; n < nvectors; n++){
               ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
               if(!inputfile) cerr << "*** ERROR 2: Lanczos could not open vectors.dat\n";
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
  dmtk::Vector<double> d(100), e(100), xc(100);
  dmtk::Matrix<double> z(100,100);
  uint i, j, col, iter = 0;
  double eini, e0;
  int control_max = maxiter;
  maxiter = 0;

  ofstream outputfile(vectors_file,std::ios::out|std::ios::binary);
  if(!outputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
       
 
  e0 = 10000.0f;
  x1 = x2 = aux = T(0);

  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    if(use_seed){
      x1 = seed;
    } else { 
      srand(9872121); // arbitrary seed
      for(i = 0; i < x1.size(); i++){
         x1(i) = random() * 1.e-10 - random() * 1.e-10; 
//         x1(i) = 1.;
      }
    } 

    iter = 0;

    b(0) = sqrt(real(product(x1,x1)));
    x1 = x1 / T(b(0));
    x2 = T(0);
    b(0) = 1.;

    uint nmax = std::min(99, (int)gs.size()-1);
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
      a(iter) = real(product(x1, x2));
      x1 = x1 - x2*T(a(iter));
      b(iter) = sqrt(real(product(x1, x1)));
   
//      cout << "Iter=" << iter << " " << a(iter) << " " << b(iter) << endl;

      if(calc_gs) x2.write(outputfile);

      if(maxiter > 0) {  // We are building the ground state;

//        gs = gs + xc(iter) * x2;
//        if(iter == maxiter) return; 

      } else{  // we are calculating energy

        if(iter >= 2){

          d(Range(1,iter)) = a(Range(1,iter));
          e(Range(2,iter+1)) = b(Range(1,iter));
          // call tqli without eigenvectors
          tqli(d, e, iter, z, false);

  
          e0 = 10000.0f;
          for(j = 1; j <= iter; j++){
            if(d(j) < e0) {
              e0 = d(j);
              col = j;
            }
          }

          cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << endl;
          if((force_maxiter && iter == control_max) || (!force_maxiter && (iter == gs.size()-1 || iter == 99 || fabs(b(iter)) < tol || fabs(eini-e0) <= tol || iter == control_max)))
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

             outputfile.flush();
             outputfile.close();

//             break; // Exit Lanczos iteration. Re-start for ground state.

         
             ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
             if(!inputfile) cerr << "*** ERROR 2: Lanczos could not open vectors.dat\n";
             gs = T(0);
             ener = d(col);
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

#undef SWAP
#undef SIGN

} // namespace dmtk

#endif // __DMTK_ALGO_H__
