#ifndef BLAS_CONNECTOR_H
#define BLAS_CONNECTOR_H

#include <complex>
using namespace std;

extern "C"
{
	// level 1: std::vector-std::vector operations	

	// Peize Lin add ?scal 2016-08-04, to compute x=a*x
	void sscal_(const int *N, const float *alpha, float *X, const int *incX);
	void dscal_(const int *N, const double *alpha, double *X, const int *incX);
	void cscal_(const int *N, const std::complex<float> *alpha, std::complex<float> *X, const int *incX);
	void zscal_(const int *N, const std::complex<double> *alpha, std::complex<double> *X, const int *incX);

	// Peize Lin add ?axpy 2016-08-04, to compute y=a*x+y
	void saxpy_(const int *N, const float *alpha, const float *X, const int *incX, float *Y, const int *incY);
	void daxpy_(const int *N, const double *alpha, const double *X, const int *incX, double *Y, const int *incY);
	void caxpy_(const int *N, const std::complex<float> *alpha, const std::complex<float> *X, const int *incX, std::complex<float> *Y, const int *incY);
	void zaxpy_(const int *N, const std::complex<double> *alpha, const std::complex<double> *X, const int *incX, std::complex<double> *Y, const int *incY);
	
	void dcopy_(long const *n, const double *a, int const *incx, double *b, int const *incy);
	void zcopy_(long const *n, const std::complex<double> *a, int const *incx, std::complex<double> *b, int const *incy); 

	void zdotc_(std::complex<double> *result, const int *n, const std::complex<double> *zx, 
		const int *incx, const std::complex<double> *zy, const int *incy);

	// Peize Lin add ?dot 2017-10-27, to compute d=x*y
	float sdot_(const int *N, const float *X, const int *incX, const float *Y, const int *incY);
	double ddot_(const int *N, const double *X, const int *incX, const double *Y, const int *incY);

	// level 2: matrix-std::vector operations
	void dgemv_(const char *transa, const int *m, const int *n, const double *alpha,  const double *a,  
		const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
		
	void zgemv_(const char *trans, const int *m, const int *n, const std::complex<double> *alpha,
			const std::complex<double> *a, const int *lda, const std::complex<double> *x, const int *incx,
			const std::complex<double> *beta, std::complex<double> *y, const int *incy);

	void dsymv_(const char *uplo, const int *n, 
		const double *alpha, const double *a, const int *lda, 
		const double *x, const int *incx, 
		const double *beta, double *y, const int *incy);			
			
	// level 3: matrix-matrix operations
	
	// Peize Lin add ?gemm 2017-10-27, to compute C = a * A.? * B.? + b * C 
	void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const float *alpha, const float *a, const int *lda, const float *b, const int *ldb, 
		const float *beta, float *c, const int *ldc);
	void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, 
		const double *beta, double *c, const int *ldc);
	void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, const std::complex<double> *b, const int *ldb, 
		const std::complex<double> *beta, std::complex<double> *c, const int *ldc);

    void dzgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
        const std::complex<double> *alpha, const double *a, const int *lda, const std::complex<double> *b,
        const int *ldb, const std::complex<double> *beta, std::complex<double> *c, const int *ldc);

	void dsymm_(const char *side, const char *uplo, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);

	void zhemm_(char *side, char *uplo, int *m, int *n,std::complex<double> *alpha,
		std::complex<double> *a,  int *lda,  std::complex<double> *b, int *ldb, std::complex<double> *beta, std::complex<double> *c, int *ldc);

	void dtrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
		double* alpha, double* a, int *lda, double*b, int *ldb);

	void ztrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
	std::complex<double>* alpha, std::complex<double>* a, int *lda, std::complex<double>*b, int *ldb);

	// computes all eigenvalues of a symmetric tridiagonal matrix
	// using the Pal-Walker-Kahan variant of the QL or QR algorithm.
	void dsterf_(int *n, double *d, double *e, int *info);

	void dstein_(int *n, double* d, double *e, int *m, double *w,
		int* block, int* isplit, double* z, int *lda, double *work,
		int* iwork, int* ifail, int *info);
 	void zstein_(int *n, double* d, double *e, int *m, double *w,
        int* block, int* isplit, std::complex<double>* z, int *lda, double *work,
        int* iwork, int* ifail, int *info);

	// computes the Cholesky factorization of a real symmetric
	// positive definite matrix A.
	void dpotf2_(char *uplo, int *n, double *a, int *lda, int *info);
	void zpotf2_(char *uplo,int *n,std::complex<double> *a, int *lda, int *info);

	void dsygs2_(int *itype, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, int *info);
	void zhegs2_(int *itype, char *uplo, int *n, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, int *info);

	void dlacpy_(char *uplo, int *m, int *n, double* a, int *lda, double *b, int *ldb);
	void zlacpy_(char *uplo, int *m, int *n, std::complex<double>* a, int *lda, std::complex<double> *b, int *ldb);

	void dlarfg_(int *n, double *alpha, double *x, int *incx, double *tau);
	void zlarfg_(int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx, std::complex<double> *tau);

	void dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
	void zgerc_(int *m, int *n, std::complex<double> *alpha,std::complex<double> *x, int *incx, std::complex<double> *y, int *incy,std::complex<double> *a, int *lda);

};

#endif
