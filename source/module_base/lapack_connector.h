#ifndef LAPACKCONNECTOR_HPP
#define LAPACKCONNECTOR_HPP

#include <new>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include "matrix.h"
#include "complexmatrix.h"
#include "global_function.h"

//Naming convention of lapack subroutines : ammxxx, where
//"a" specifies the data type:
//  - d stands for double
//  - z stands for complex double
//"mm" specifies the type of matrix, for example:
//  - he stands for hermitian
//  - sy stands for symmetric
//"xxx" specifies the type of problem, for example:
//  - gv stands for generalized eigenvalue

extern "C"
{
    int ilaenv_(int* ispec,const char* name,const char* opts,
    const int* n1,const int* n2,const int* n3,const int* n4);
    // solve the generalized eigenproblem Ax=eBx, where A is Hermitian and complex couble
    // zhegv_ returns all eigenvalues while zhegvx_ returns selected ones
    void zhegv_(const int* itype,const char* jobz,const char* uplo,const int* n,
                std::complex<double>* a,const int* lda,std::complex<double>* b,const int* ldb,
                double* w,std::complex<double>* work,int* lwork,double* rwork,int* info);
    void zhegvx_(const int* itype,const char* jobz,const char* range,const char* uplo,
                 const int* n,std::complex<double> *a,const int* lda,std::complex<double> *b,
                 const int* ldb,const double* vl,const double* vu,const int* il,
                 const int* iu,const double* abstol,const int* m,double* w,
                 std::complex<double> *z,const int *ldz,std::complex<double> *work,const int* lwork,
                 double* rwork,int* iwork,int* ifail,int* info);
	// solve the generalized eigenproblem Ax=eBx, where A is Symmetric and real couble
    // dsygv_ returns all eigenvalues while dsygvx_ returns selected ones
	void dsygv_(const int* itype, const char* jobz,const char* uplo, const int* n,
				double* a,const int* lda,double* b,const int* ldb,
	 			double* w,double* work,int* lwork,int* info);
	void dsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
                const int* n, double* A, const int* lda, double* B, const int* ldb,
                const double* vl, const double* vu, const int* il, const int* iu,
                const double* abstol, int* m, double* w, double* Z, const int* ldz, 
                double* work, int* lwork, int*iwork, int* ifail, int* info);
    // solve the eigenproblem Ax=ex, where A is Symmetric and real double
	void dsyev_(const char* jobz,const char* uplo,const int* n,double *a,
                const int* lda,double* w,double* work,const int* lwork, int* info);
    // solve the eigenproblem Ax=ex, where A is Hermitian and complex couble
    void zheev_(const char* jobz,const char* uplo,const int* n,std::complex<double> *a,
                const int* lda,double* w,std::complex<double >* work,const int* lwork,
                double* rwork,int* info);
    // dsytrf_ computes the Bunch-Kaufman factorization of a double precision
    // symmetric matrix, while dsytri takes its output to perform martrix inversion
    void dsytrf_(const char* uplo, const int* n, double * a, const int* lda,
                 int *ipiv,double *work, int* lwork ,int *info);
    void dsytri_(const char* uplo,const int* n,double *a, const int *lda,
                 int *ipiv, double * work,int *info);
    // Peize Lin add dsptrf and dsptri 2016-06-21, to compute inverse real symmetry indefinit matrix.
    // dpotrf computes the Cholesky factorization of a real symmetric positive definite matrix
    // while dpotri taks its output to perform matrix inversion
    void dpotrf_(const char* uplo,const int* n, double* A, const int* lda, int *info);
    void dpotri_(const char* uplo,const int* n, double* A, const int* lda, int *info);

    // zgetrf computes the LU factorization of a general matrix
    // while zgetri takes its output to perform matrix inversion
    void zgetrf_(const int* m, const int *n, const std::complex<double> *A, const int *lda, int *ipiv, const int* info);
    void zgetri_(const int* n, std::complex<double> *A, const int *lda, int *ipiv, std::complex<double> *work, int *lwork, const int *info);

    // if trans=='N':	C = alpha * A * A.H + beta * C
	// if trans=='C':	C = alpha * A.H * A + beta * C
	void zherk_(const char *uplo, const char *trans, const int *n, const int *k, 
		const double *alpha, const std::complex<double> *A, const int *lda, 
		const double *beta, std::complex<double> *C, const int *ldc);

	// computes all eigenvalues of a symmetric tridiagonal matrix
	// using the Pal-Walker-Kahan variant of the QL or QR algorithm.
	void dsterf_(int *n, double *d, double *e, int *info);
    // computes the eigenvectors of a real symmetric tridiagonal
    // matrix T corresponding to specified eigenvalues
	void dstein_(int *n, double* d, double *e, int *m, double *w,
		int* block, int* isplit, double* z, int *lda, double *work,
		int* iwork, int* ifail, int *info);
    // computes the eigenvectors of a complex symmetric tridiagonal
    // matrix T corresponding to specified eigenvalues
 	void zstein_(int *n, double* d, double *e, int *m, double *w,
        int* block, int* isplit, std::complex<double>* z, int *lda, double *work,
        int* iwork, int* ifail, int *info);

	// computes the Cholesky factorization of a symmetric
	// positive definite matrix A.
	void dpotf2_(char *uplo, int *n, double *a, int *lda, int *info);
	void zpotf2_(char *uplo,int *n,std::complex<double> *a, int *lda, int *info);

    // reduces a symmetric definite generalized eigenproblem to standard form
    // using the factorization results obtained from spotrf
	void dsygs2_(int *itype, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, int *info);
	void zhegs2_(int *itype, char *uplo, int *n, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, int *info);

    // copies a into b
	void dlacpy_(char *uplo, int *m, int *n, double* a, int *lda, double *b, int *ldb);
	void zlacpy_(char *uplo, int *m, int *n, std::complex<double>* a, int *lda, std::complex<double> *b, int *ldb);

    // generates a real elementary reflector H of order n, such that
    //   H * ( alpha ) = ( beta ),   H is unitary.
    //       (   x   )   (   0  )
	void dlarfg_(int *n, double *alpha, double *x, int *incx, double *tau);
	void zlarfg_(int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx, std::complex<double> *tau);

}

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
class LapackConnector
{
private:
    // Transpose the std::complex matrix to the fortran-form real-std::complex array.
    static inline
    std::complex<double>* transpose(const ModuleBase::ComplexMatrix& a, const int n, const int lda)
    {
        std::complex<double>* aux = new std::complex<double>[lda*n];
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < lda; ++j)
            {
                aux[i*lda+j] = a(j,i);		// aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
            }
        }
        return aux;
    }

    // Transpose the fortran-form real-std::complex array to the std::complex matrix.
    static inline
    void transpose(const std::complex<double>* aux, ModuleBase::ComplexMatrix& a, const int n, const int lda)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < lda; ++j)
            {
                a(j, i) = aux[i*lda+j];		// aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
            }
        }
    }
	
	// Peize Lin add 2015-12-27	
	static inline
	char change_uplo(const char &uplo)
	{
		switch(uplo)
		{
			case 'U': return 'L';
			case 'L': return 'U';
			default: throw std::invalid_argument("uplo must be 'U' or 'L'");
		}		
	}
	
	// Peize Lin add 2019-04-14
	static inline
	char change_trans_NC(const char &trans)
	{
		switch(trans)
		{
			case 'N': return 'C';
			case 'C': return 'N';
			default: throw std::invalid_argument("trans must be 'N' or 'C'");
		}
	}
	
public:

    static inline
    int ilaenv( int ispec, const char *name,const char *opts,const int n1,const int n2,
                const int n3,const int n4)
    {
        const int nb = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);
        return nb;
    }

    // wrap function of fortran lapack routine zhegv.
    static inline
    void zhegv(	const int itype,const char jobz,const char uplo,const int n,ModuleBase::ComplexMatrix& a,
                const int lda,ModuleBase::ComplexMatrix& b,const int ldb,double* w,std::complex<double>* work,
                int lwork,double* rwork,int info)
    {	// Transpose the std::complex matrix to the fortran-form real-std::complex array.
        std::complex<double>* aux = LapackConnector::transpose(a, n, lda);
        std::complex<double>* bux = LapackConnector::transpose(b, n, ldb);

        // call the fortran routine
        zhegv_(&itype, &jobz, &uplo, &n, aux, &lda, bux, &ldb, w, work, &lwork, rwork, &info);
        // Transpose the fortran-form real-std::complex array to the std::complex matrix.
        LapackConnector::transpose(aux, a, n, lda);
        LapackConnector::transpose(bux, b, n, ldb);
        // free the memory.
        delete[] aux;
        delete[] bux;
    }

    // wrap function of fortran lapack routine zheev.
    static inline
    void zhegvx( const int itype, const char jobz, const char range, const char uplo,
                 const int n, const ModuleBase::ComplexMatrix& a, const int lda, const ModuleBase::ComplexMatrix& b,
                 const int ldb, const double vl, const double vu, const int il, const int iu,
                 const double abstol, const int m, double* w, ModuleBase::ComplexMatrix& z, const int ldz,
                 std::complex<double>* work, const int lwork, double* rwork, int* iwork,
                 int* ifail, int info)
    {
        // Transpose the std::complex matrix to the fortran-form real-std::complex array.
        std::complex<double>* aux = LapackConnector::transpose(a, n, lda);
        std::complex<double>* bux = LapackConnector::transpose(b, n, ldb);
        std::complex<double>* zux = new std::complex<double>[n*iu];// mohan modify 2009-08-02
        //for(int i=0; i<n*iu; i++) zux[i] = std::complex<double>(0.0,0.0);

        // call the fortran routine
        zhegvx_(&itype, &jobz, &range, &uplo, &n, aux, &lda, bux, &ldb, &vl,
                &vu, &il,&iu, &abstol, &m, w, zux, &ldz, work, &lwork, rwork, iwork, ifail, &info);

        // Transpose the fortran-form real-std::complex array to the std::complex matrix
        LapackConnector::transpose(zux, z, iu, n);	 // mohan modify 2009-08-02

        // free the memory.
        delete[] aux;
        delete[] bux;
        delete[] zux;

    }	

	// calculate the eigenvalues and eigenfunctions of a real symmetric matrix.
    static inline
    void dsygv(	const int itype,const char jobz,const char uplo,const int n,ModuleBase::matrix& a,
                const int lda,ModuleBase::matrix& b,const int ldb,double* w,double* work,
                int lwork,int *info)
    {	// Transpose the std::complex matrix to the fortran-form real-std::complex array.
        double* aux = new double[lda*n];
        double* bux = new double[lda*n];
	    for (int i=0; i<n; i++)
        {
            for (int j=0; j<lda; j++)
            {
                aux[i*lda+j] = a(j,i);      // aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
				bux[i*lda+j] = b(j,i);
            }
        }
        // call the fortran routine
        dsygv_(&itype, &jobz, &uplo, &n, aux, &lda, bux, &ldb, w, work, &lwork, info);
		for (int i=0; i<n; i++)
		{
			for(int j=0; j<lda; ++j)
			{
				a(j,i)=aux[i*lda+j];
				b(j,i)=bux[i*lda+j];
			}
		}
 
        // free the memory.
        delete[] aux;
        delete[] bux;
    }

	// wrap function of fortran lapack routine dsyev.
    static inline
	void dsyev( const char jobz, const char uplo, ModuleBase::matrix &a, double* w, int info )		// Peize Lin update 2017-10-17
	{
		assert(a.nr==a.nc);
		const char uplo_changed = change_uplo(uplo);
		
		double work_tmp;
		constexpr int minus_one = -1;
		dsyev_(&jobz, &uplo_changed, &a.nr, a.c, &a.nr, w, &work_tmp, &minus_one, &info);		// get best lwork
		
		const int lwork = work_tmp;
		std::vector<double> work(std::max(1,lwork));
		dsyev_(&jobz, &uplo_changed, &a.nr, a.c, &a.nr, w, ModuleBase::GlobalFunc::VECTOR_TO_PTR(work), &lwork, &info);		
	}
    
    // wrap function of fortran lapack routine zheev.
    static inline
    void zheev( const char jobz,
                const char uplo,
                const int n,
                ModuleBase::ComplexMatrix& a,
                const int lda,
                double* w,
                std::complex< double >* work,
                const int lwork,
                double* rwork,
                int *info	)
    {	// Transpose the std::complex matrix to the fortran-form real-std::complex array.
        std::complex<double> *aux = LapackConnector::transpose(a, n, lda);
        // call the fortran routine
        zheev_(&jobz, &uplo, &n, aux, &lda, w, work, &lwork, rwork, info);
        // Transpose the fortran-form real-std::complex array to the std::complex matrix.
        LapackConnector::transpose(aux, a, n, lda);
        // free the memory.
        delete[] aux;
    }

    static inline
    void zgetrf(int m, int n, ModuleBase::ComplexMatrix &a, const int lda, int *ipiv, int *info)
    {
        std::complex<double> *aux = LapackConnector::transpose(a, n, lda);
        zgetrf_( &m, &n, aux, &lda, ipiv, info);
        LapackConnector::transpose(aux, a, n, lda);
        delete[] aux;
                return;
    }
    static inline
    void zgetri(int n, ModuleBase::ComplexMatrix &a,  int lda, int *ipiv, std::complex<double> * work, int lwork, int *info)
    {
        std::complex<double> *aux = LapackConnector::transpose(a, n, lda);
        zgetri_( &n, aux, &lda, ipiv, work, &lwork, info);
        LapackConnector::transpose(aux, a, n, lda);
        delete[] aux;
        return;
    }

	// Peize Lin add 2016-07-09
	static inline
	void dpotrf( char uplo, const int n, ModuleBase::matrix &a, const int lda, int *info )
	{
		const char uplo_changed = change_uplo(uplo);
		dpotrf_( &uplo_changed, &n, a.c, &lda, info );
	}	
	
	// Peize Lin add 2016-07-09
	static inline
	void dpotri( char uplo, const int n, ModuleBase::matrix &a, const int lda, int *info )
	{
		const char uplo_changed = change_uplo(uplo);
		dpotri_( &uplo_changed, &n, a.c, &lda, info);		
	}	
	
	// Peize Lin add 2019-04-14
	// if trans=='N':	C = a * A * A.H + b * C
	// if trans=='C':	C = a * A.H * A + b * C
	static inline
	void zherk(const char uplo, const char trans, const int n, const int k, 
		const double alpha, const std::complex<double> *A, const int lda, 
		const double beta, std::complex<double> *C, const int ldc)
	{
		const char uplo_changed = change_uplo(uplo);
		const char trans_changed = change_trans_NC(trans);
		zherk_(&uplo_changed, &trans_changed, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
	}
	
};
#endif  // LAPACKCONNECTOR_HPP
