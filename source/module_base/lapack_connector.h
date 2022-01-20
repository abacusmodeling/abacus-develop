#ifndef LAPACKCONNECTOR_HPP
#define LAPACKCONNECTOR_HPP

#include <new>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include "matrix.h"
#include "complexmatrix.h"
#include "global_function.h"

extern "C"
{
    int ilaenv_(int* ispec,const char* name,const char* opts,
    const int* n1,const int* n2,const int* n3,const int* n4);
    void zhegv_(const int* itype,const char* jobz,const char* uplo,const int* n,
                std::complex<double>* a,const int* lda,std::complex<double>* b,const int* ldb,
                double* w,std::complex<double>* work,int* lwork,double* rwork,int* info);
	// mohan add dsygv 2010-03-20, to compute the eigenfunctions and eigenvalues of a real symmetry matrix.
	void dsygv_(const int* itype, const char* jobz,const char* uplo, const int* n,
				double* a,const int* lda,double* b,const int* ldb,
	 			double* w,double* work,int* lwork,int* info);
    // Peize Lin add dsygvx 2020.11.11, to compute the selected eigenvalues of a generalized symmetric/Hermitian matrix.
	void dsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
                const int* n, double* A, const int* lda, double* B, const int* ldb,
                const double* vl, const double* vu, const int* il, const int* iu,
                const double* abstol, int* m, double* w, double* Z, const int* ldz, 
                double* work, int* lwork, int*iwork, int* ifail, int* info);
	void dsyev_(const char* jobz,const char* uplo,const int* n,double *a,
                const int* lda,double* w,double* work,const int* lwork, int* info);
    void zheev_(const char* jobz,const char* uplo,const int* n,std::complex<double> *a,
                const int* lda,double* w,std::complex<double >* work,const int* lwork,
                double* rwork,int* info);
    void dsytrf_(const char* uplo, const int* n, double * a, const int* lda,
                 int *ipiv,double *work, int* lwork ,int *info);
    void dsytri_(const char* uplo,const int* n,double *a, const int *lda,
                 int *ipiv, double * work,int *info);
    void zhegvx_(const int* itype,const char* jobz,const char* range,const char* uplo,
                 const int* n,std::complex<double> *a,const int* lda,std::complex<double> *b,
                 const int* ldb,const double* vl,const double* vu,const int* il,
                 const int* iu,const double* abstol,const int* m,double* w,
                 std::complex<double> *z,const int *ldz,std::complex<double> *work,const int* lwork,
                 double* rwork,int* iwork,int* ifail,int* info);
    // Peize Lin add dsptrf and dsptri 2016-06-21, to compute inverse real symmetry indefinit matrix.
    void dpotrf_(const char* uplo,const int* n, double* A, const int* lda, int *info);
    void dpotri_(const char* uplo,const int* n, double* A, const int* lda, int *info);

    void zgetrf_(const int* m, const int *n, const std::complex<double> *A, const int *lda, int *ipiv, const int* info);
    void zgetri_(const int* n, std::complex<double> *A, const int *lda, int *ipiv, std::complex<double> *work, int *lwork, const int *info);
	void zherk_(const char *uplo, const char *trans, const int *n, const int *k, 
		const double *alpha, const std::complex<double> *A, const int *lda, 
		const double *beta, std::complex<double> *C, const int *ldc);
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
	void dsyev(const char jobz,const char uplo,const int n,double *a,
			const int lda,double *w,double *work,const int lwork, int info)
	{
		const char uplo_changed = change_uplo(uplo);
		dsyev_(&jobz, &uplo_changed, &n, a, &lda, w, work, &lwork, &info);
	}
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
