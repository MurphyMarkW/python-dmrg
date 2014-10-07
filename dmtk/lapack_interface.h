#ifndef __DMTK_LAPACK_INTERFACE_H__
#define __DMTK_LAPACK_INTERFACE_H__

#include <complex>

extern "C" {

#if defined(MKL_ILP64) || defined(__LP64__) || defined(__APPLE__)
typedef long DMTK_int;
typedef double DMTK_double;
#else
typedef int DMTK_int;
typedef double DMTK_double;
#endif

void dgesvd_(
        const char &jobu,               // (input)
        const char &jobvt,              // (input)
        const DMTK_int &m,                   // (input)
        const DMTK_int &n,                   // (input)
        double *a,                      // a[n][lda] (input/output)
        const DMTK_int &lda,                 // (input)
        double *s,                      // s[min(m,n)] (output)
        double *u,                      // u[ucol][ldu] (output)
        const DMTK_int &ldu,                 // (input)
        double *vt,                     // vt[n][ldvt] (output)
        const DMTK_int &ldvt,                // (input)
        double *work,                   // work[lwork] (workspace/output)
        const DMTK_int &lwork,               // (input)
        DMTK_int &info                       // (output)
        );

void zheev_(
	const char &jobz,		// (input)
	const char &uplo,		// (input)
	const DMTK_int &n,			// (input)
	std::complex<double> *a,	// a[n][lda] (input/output)
	const DMTK_int &lda,			// (input)
	double *w,			// w[n] (output)
	std::complex<double> *work,	// work[lwork] (workspace/output)
	const DMTK_int &lwork,		// (input)
	double *rwork,			// rwork[max(1, 3*n-2)] (workspace)
	DMTK_int &info			// (output)
	);

void ssyev_(
	const char &jobz,		// (input)
	const char &uplo,		// (input)
	const DMTK_int &n,			// (input)
	float *a,			// a[n][lda] (input/output)
	const DMTK_int &lda,			// (input)
	float *w,			// w[n] (output)
	float *work,			// work[lwork] (workspace/output)
	const DMTK_int &lwork,		// (input)
	DMTK_int &info			// (output)
	);

void dsyev_(
	const char &jobz,		// (input)
	const char &uplo,		// (input)
	const DMTK_int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const DMTK_int &lda,			// (input)
	double *w,			// w[n] (output)
	double *work,			// work[lwork] (workspace/output)
	const DMTK_int &lwork,		// (input)
	DMTK_int &info			// (output)
	);

void zhpev_(
	const char &jobz,		// (input)
	const char &uplo,		// (input)
	const DMTK_int &n,			// (input)
	std::complex<double> *ap,	// ap[n*(n+1)/2] (input/output)
	double *w,			// w[n] (output)
	std::complex<double> *z,	// z[n][ldz] (output)
	const DMTK_int &ldz,			// (input)
	std::complex<double> *work,	// work[max(1, 2*n-1)] (workspace)
	double *rwork,			// rwork[max(1, 3*n-2)] (workspace)
	DMTK_int &info			// (output)
	);

void dstevd_(
	const char &jobz,		// (input)
	const DMTK_int &n,			// (input)
	double *d,			// d[n] (input/output)
	double *e,			// e[n] (input/output)
	double *z,			// z[n][ldz] (output)
	const DMTK_int &ldz,			// (input)
	double *work,			// work[?] (workspace/output)
	const DMTK_int &lwork,		// (input)
	DMTK_int *iwork,			// iwork[liwork] (workspace/output)
	const DMTK_int &liwork,		// (input)
	DMTK_int &info			// (output)
	);

void dgesv_(const DMTK_int & n,
            const DMTK_int & nrhs,
            double *da,
            const DMTK_int & lda,
            DMTK_int *ipivot,
            double *db,
            const DMTK_int & ldb,
            DMTK_int & info);

void sgesv_(const DMTK_int & n,
            const DMTK_int & nrhs,
            float *sa,
            const DMTK_int & lda,
            DMTK_int *ipivot,
            float *sb,
            const DMTK_int & ldb,
            DMTK_int & info);

void zgesv_(const DMTK_int & n,
            const DMTK_int & nrhs,
            std::complex<double> *za,
            const DMTK_int & lda,
            DMTK_int *ipivot,
            std::complex<double> *zb,
            const DMTK_int & ldb,
            DMTK_int & info);

void cgesv_(const DMTK_int & n,
            const DMTK_int & nrhs,
            std::complex<float> *ca,
            const DMTK_int & lda,
            DMTK_int *ipivot,
            std::complex<float> *cb,
            const DMTK_int & ldb,
            DMTK_int & info);

void dgeevx_(
	const char &balanc,		// (input)
	const char &jobvl,		// (input)
	const char &jobvr,		// (input)
	const char &sense,		// (input)
	const DMTK_int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const DMTK_int &lda,			// (input)
	double *wr,			// wr[n] (output)
	double *wi,			// wi[n] (output)
	double *vl,			// vl[n][ldvl] (output)
	const DMTK_int &ldvl,		// (input)
	double *vr,			// vr[n][ldvr] (output)
	const DMTK_int &ldvr,		// (input)
	DMTK_int &ilo,			// (output)
	DMTK_int &ihi,			// (output)
	double *scale,			// scale[n] (output)
	double &abnrm,			// (output)
	double *rconde,			// rconde[n] (output)
	double *rcondv,			// rcondv[n] (output)
	double *work,			// work[lwork] (workspace/output)
	const DMTK_int &lwork,		// (input)
	DMTK_int *iwork,			// iwork[2*n-2] (workspace)
	DMTK_int &info			// (output)
	);

double sdot_(
            const DMTK_int&,
            const float*, const DMTK_int&,
            const float*, const DMTK_int&);

double ddot_(
            const DMTK_int&,
            const double*, const DMTK_int&,
            const double*, const DMTK_int&);

std::complex<double> zdotc_(
            const DMTK_int&,
            const std::complex<double>*, const DMTK_int&,
            const std::complex<double>*, const DMTK_int&);

void dcopy_(const DMTK_int &,
            const double *, const DMTK_int &,
            double *, const DMTK_int &);

void zcopy_(const DMTK_int &,
            const std::complex<double> *, const DMTK_int &,
            std::complex<double> *, const DMTK_int &);

void sgemv_(const char&,
            const DMTK_int&, const DMTK_int&,
            const float&,
            const float*, const DMTK_int&,
            const float*, const DMTK_int&,
            const float&,
            float*, const DMTK_int&);

void dgemv_(const char&,
            const DMTK_int&, const DMTK_int&,
            const double&,
            const double*, const DMTK_int&,
            const double*, const DMTK_int&,
            const double&,
            double*, const DMTK_int&);

void zgemv_(const char&,
            const DMTK_int&, const DMTK_int&,
            const std::complex<double>&,
            const std::complex<double>*, const DMTK_int&,
            const std::complex<double>*, const DMTK_int&,
            const std::complex<double>&,
            std::complex<double>*, const DMTK_int&);

void dgemm_(const char&, const char&,
            const DMTK_int&, const DMTK_int&, const DMTK_int&,
            const double&,
            const double*, const DMTK_int&,
            const double*, const DMTK_int& ,
            const double&,
            double *, const DMTK_int&);

void zgemm_(const char&, const char&,
            const DMTK_int&, const DMTK_int&, const DMTK_int&,
            const std::complex<double>&,
            const std::complex<double>*, const DMTK_int&,
            const std::complex<double>*, const DMTK_int&,
            const std::complex<double>&,
            std::complex<double>*, const DMTK_int&);

void sgemm_(const char&, const char&,
            const DMTK_int&, const int&, const DMTK_int&,
            const float&,
            const float*, const DMTK_int&,
            const float*, const DMTK_int& ,
            const float&,
            float *, const DMTK_int&);

void daxpy_(const DMTK_int &n, const double& da,
            const double* dx, const DMTK_int& incx,
            const double* dy, const DMTK_int& incy);

void saxpy_(const int &n, const float& da,
            const float* dx, const DMTK_int& incx,
            const float* dy, const DMTK_int& incy);

}

#endif // __DMTK_LAPACK_INTERFACE_H__
