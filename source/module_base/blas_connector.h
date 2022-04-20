#ifndef BLAS_CONNECTOR_H
#define BLAS_CONNECTOR_H

#include <complex>

extern "C"
{
	// level 1: std::vector-std::vector operations, O(n) data and O(n) work.

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

	//reason for passing results as argument instead of returning it:
	//see https://www.numbercrunch.de/blog/2014/07/lost-in-translation/
	void zdotc_(std::complex<double> *result, const int *n, const std::complex<double> *zx, 
		const int *incx, const std::complex<double> *zy, const int *incy);
	// Peize Lin add ?dot 2017-10-27, to compute d=x*y
	float sdot_(const int *N, const float *X, const int *incX, const float *Y, const int *incY);
	double ddot_(const int *N, const double *X, const int *incX, const double *Y, const int *incY);

	// Peize Lin add ?nrm2 2018-06-12, to compute out = ||x||_2 = \sqrt{ \sum_i x_i**2 }
	float snrm2_( const int *n, const float *X, const int *incX );
	double dnrm2_( const int *n, const double *X, const int *incX );
	double dznrm2_( const int *n, const std::complex<double> *X, const int *incX );

	// level 2: matrix-std::vector operations, O(n^2) data and O(n^2) work.
	void dgemv_(const char *transa, const int *m, const int *n, const double *alpha,  const double *a,  
		const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
		
	void zgemv_(const char *trans, const int *m, const int *n, const std::complex<double> *alpha,
			const std::complex<double> *a, const int *lda, const std::complex<double> *x, const int *incx,
			const std::complex<double> *beta, std::complex<double> *y, const int *incy);

	void dsymv_(const char *uplo, const int *n, 
		const double *alpha, const double *a, const int *lda, 
		const double *x, const int *incx, 
		const double *beta, double *y, const int *incy);			

    // A := alpha x * y.T + A
	void dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
	void zgerc_(int *m, int *n, std::complex<double> *alpha,std::complex<double> *x, int *incx, std::complex<double> *y, int *incy,std::complex<double> *a, int *lda);

	// level 3: matrix-matrix operations, O(n^2) data and O(n^3) work.
	
	// Peize Lin add ?gemm 2017-10-27, to compute C = a * A.? * B.? + b * C
	// A is general
	void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const float *alpha, const float *a, const int *lda, const float *b, const int *ldb, 
		const float *beta, float *c, const int *ldc);
	void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, 
		const double *beta, double *c, const int *ldc);
	void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, const std::complex<double> *b, const int *ldb, 
		const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
	
	//a is symmetric
	void dsymm_(const char *side, const char *uplo, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
	//a is hermitian
	void zhemm_(char *side, char *uplo, int *m, int *n,std::complex<double> *alpha,
		std::complex<double> *a,  int *lda,  std::complex<double> *b, int *ldb, std::complex<double> *beta, std::complex<double> *c, int *ldc);

	//solving triangular matrix with multiple right hand sides
	void dtrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
		double* alpha, double* a, int *lda, double*b, int *ldb);
	void ztrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
	std::complex<double>* alpha, std::complex<double>* a, int *lda, std::complex<double>*b, int *ldb);

};

// Class BlasConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	BlasConnector::functionname(parameter list).
class BlasConnector
{
public:

	// Peize Lin add 2016-08-04
	// y=a*x+y
	static inline 
	void axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY)
	{
		saxpy_(&n, &alpha, X, &incX, Y, &incY);
	}	
	static inline 
	void axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY)
	{
		daxpy_(&n, &alpha, X, &incX, Y, &incY);
	}	
	static inline 
	void axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY)
	{
		caxpy_(&n, &alpha, X, &incX, Y, &incY);
	}	
	static inline 
	void axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY)
	{
		zaxpy_(&n, &alpha, X, &incX, Y, &incY);
	}	
	
	// Peize Lin add 2016-08-04
	// x=a*x
	static inline 
	void scal( const int n,  const float alpha, float *X, const int incX)
	{
		sscal_(&n, &alpha, X, &incX);
	}	
	static inline 
	void scal( const int n, const double alpha, double *X, const int incX)
	{
		dscal_(&n, &alpha, X, &incX);
	}	
	static inline 
	void scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX)
	{
		cscal_(&n, &alpha, X, &incX);
	}	
	static inline 
	void scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX)
	{
		zscal_(&n, &alpha, X, &incX);
	}	
	
	// Peize Lin add 2017-10-27
	// d=x*y
	static inline
	float dot( const int n, const float *X, const int incX, const float *Y, const int incY)
	{
		return sdot_(&n, X, &incX, Y, &incY);
	}
	static inline
	double dot( const int n, const double *X, const int incX, const double *Y, const int incY)
	{
		return ddot_(&n, X, &incX, Y, &incY);
	}

	// Peize Lin add 2017-10-27, fix bug trans 2019-01-17
	// C = a * A.? * B.? + b * C 
	static inline
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const float alpha, const float *a, const int lda, const float *b, const int ldb, 
		const float beta, float *c, const int ldc)
	{
		sgemm_(&transb, &transa, &n, &m, &k,
			&alpha, b, &ldb, a, &lda, 
			&beta, c, &ldc);
	}
	static inline
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const double alpha, const double *a, const int lda, const double *b, const int ldb, 
		const double beta, double *c, const int ldc)
	{
		dgemm_(&transb, &transa, &n, &m, &k,
			&alpha, b, &ldb, a, &lda, 
			&beta, c, &ldc);
	}
	static inline
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb, 
		const std::complex<double> beta, std::complex<double> *c, const int ldc)
	{
		zgemm_(&transb, &transa, &n, &m, &k,
			&alpha, b, &ldb, a, &lda, 
			&beta, c, &ldc);
	}

	// Peize Lin add 2018-06-12
	// out = ||x||_2
	static inline
	float nrm2( const int n, const float *X, const int incX )
	{
		return snrm2_( &n, X, &incX );
	}
	static inline
	double nrm2( const int n, const double *X, const int incX )
	{
		return dnrm2_( &n, X, &incX );
	}
	static inline
	double nrm2( const int n, const std::complex<double> *X, const int incX )
	{
		return dznrm2_( &n, X, &incX );
	}

	// copies a into b	
	static inline
	void copy(const long n, const double *a, const int incx, double *b, const int incy)
	{
		dcopy_(&n, a, &incx, b, &incy);
	}
	static inline
	void copy(const long n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy)
	{
		zcopy_(&n, a, &incx, b, &incy);
	}	

};

#endif
