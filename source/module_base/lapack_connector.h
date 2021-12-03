#ifndef LAPACKCONNECTOR_HPP
#define LAPACKCONNECTOR_HPP

#include <new>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include "matrix.h"
#include "complexmatrix.h"
#include "blas_connector.h"
#include "global_function.h"


extern "C"
{
    int ilaenv_(int* ispec,const char* name,const char* opts,
    const int* n1,const int* n2,const int* n3,const int* n4);
    void zgesv_(const int* n,const int* nrhs,const std::complex<double> *A,const int *lda,
                const int *ipiv,const std::complex<double> *B,const int* ldb,int* info);
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
	// mohan add sspgvx 2010-03-21, to compute the selected eigenvalues of a generalized symmetric/Hermitian 
	// definite generalized problems. 
	void sspgvx_(const int* itype, const char* jobz,const char* range,const char* uplo,
				const int* n,double *ap,double *bp,const double* vl,
				const int* vu, const int* il, const int* iu,const double* abstol,
				const int* m,double* w,double *z,const int *ldz,
				double *work,int* iwork,int* ifail,int* info);
	void dsyev_(const char* jobz,const char* uplo,const int* n,double *a,
                const int* lda,double* w,double* work,const int* lwork, int* info);
    void zheev_(const char* jobz,const char* uplo,const int* n,std::complex<double> *a,
                const int* lda,double* w,std::complex<double >* work,const int* lwork,
                double* rwork,int* info);
    void dsytrf_(const char* uplo, const int* n, double * a, const int* lda,
                 int *ipiv,double *work, int* lwork ,int *info);
    void dsytri_(const char* uplo,const int* n,double *a, const int *lda,
                 int *ipiv, double * work,int *info);
    void dgesvd_(const char* jobu,const char* jobvt,const int *m,const int *n,double *a,
                 const int *lda,double *s,double *u,const int *ldu,double *vt,
                 const int *ldvt,double *work,const int *lwork,int *info);
    void zhegvx_(const int* itype,const char* jobz,const char* range,const char* uplo,
                 const int* n,std::complex<double> *a,const int* lda,std::complex<double> *b,
                 const int* ldb,const double* vl,const double* vu,const int* il,
                 const int* iu,const double* abstol,const int* m,double* w,
                 std::complex<double> *z,const int *ldz,std::complex<double> *work,const int* lwork,
                 double* rwork,int* iwork,int* ifail,int* info);
    // mohan add dsptrf and dsptri 2010-02-03, to compute inverse real symmetry indefinit matrix.
    void spotrf_(const char* uplo,const int* n, double* A, int* lda, int *info);
    void spotri_(const char* uplo,const int* n, double* A, int* lda, int *info);
    // Peize Lin add dsptrf and dsptri 2016-06-21, to compute inverse real symmetry indefinit matrix.
    void dpotrf_(const char* uplo,const int* n, double* A, const int* lda, int *info);
    void dpotri_(const char* uplo,const int* n, double* A, const int* lda, int *info);
    // mohan add zpotrf and zpotri 2010-02-03, to compute inverse std::complex hermit matrix.
    void zpotrf_(const char* uplo,const int* n,const std::complex<double> *A,const int* lda,const int* info);
    void zpotri_(const char* uplo,const int* n,const std::complex<double> *A,const int* lda,const int* info);
    void zgetrf_(const int* m, const int *n, const std::complex<double> *A, const int *lda, int *ipiv, const int* info);
    void zgetri_(const int* n, std::complex<double> *A, const int *lda, int *ipiv, std::complex<double> *work, int *lwork, const int *info);
    void pdgetrf_(const int *m, const int *n, double *a, const int *ia, const int *ja, int *desca, int *ipiv, int *info);
	// Peize Lin add ?nrm2 2018-06-12, to compute out = ||x||_2 = \sqrt{ \sum_i x_i**2 }
	float snrm2_( const int *n, const float *X, const int *incX );
	double dnrm2_( const int *n, const double *X, const int *incX );
	double dznrm2_( const int *n, const std::complex<double> *X, const int *incX );
	// Peize Lin add zherk 2019-04-14
	// if trans=='N':	C = a * A * A.H + b * C
	// if trans=='C':	C = a * A.H * A + b * C	
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
    // Transpose the  matrix to the fortran-form real array.
    static inline
    double * transpose_matrix(const ModuleBase::matrix& a, const int n, const int lda)
    {
        double * aux = new double [lda*n];
		std::cout << " lda=" << lda << " n=" << n << std::endl;
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<lda; j++)
            {
                aux[i*lda+j] = a(j,i);      // aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
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
    
	// Transpose the fortran-form real array to the matrix.
    static inline
    void transpose_matrix(const double *aux, ModuleBase::matrix &a, const int n, const int lda)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < lda; ++j)
            {
                a(j, i) = aux[i*lda+j];     // aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
            }
        }
    }
	
	// Peize Lin add 2015-12-27
	static inline
	char change_side(const char &side)
	{
		switch(side)
		{
			case 'L': return 'R';
			case 'R': return 'L';
			default: throw std::invalid_argument("Side must be 'L' or 'R'");
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
	char change_trans_NT(const char &trans)
	{
		switch(trans)
		{
			case 'N': return 'T';
			case 'T': return 'N';
			default: throw std::invalid_argument("trans must be 'N' or 'T'");
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

    static inline
    void zgesv (const int n,const int nrhs,ModuleBase::ComplexMatrix &A,const int lda,
                const int *ipiv,std::complex<double> *B,const int ldb,int* info)
    {
        std::complex<double> *aux = LapackConnector::transpose(A, n, lda);
        zgesv_(&n, &nrhs, aux, &lda, ipiv, B, &ldb, info);
        LapackConnector::transpose(aux, A, n, lda);
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


	// calculate the selected eigenvalues.
	// mohan add 2010/03/21
	static inline
	void sspgvx(const int itype, const char jobz,const char range,const char uplo,
				const int n,const ModuleBase::matrix &ap,const ModuleBase::matrix &bp,const double vl,
				const int vu, const int il, const int iu,const double abstol,
				const int m,double* w,ModuleBase::matrix &z,const int ldz,
				double *work,int* iwork,int* ifail,int* info)
	{
	    // Transpose the std::complex matrix to the fortran-form real*16 array.
		double* aux = LapackConnector::transpose_matrix(ap, n, n);
		double* bux = LapackConnector::transpose_matrix(bp, n, n);
		double* zux = new double[n*iu];

        // call the fortran routine
        sspgvx_(&itype, &jobz, &range, &uplo, 
				&n, aux, bux, &vl,
                &vu, &il,&iu, &abstol, 
				&m, w, zux, &ldz, 
				work, iwork, ifail, info);

        // Transpose the fortran-form real*16 array to the matrix
        LapackConnector::transpose_matrix(zux, z, iu, n);

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
	
    // wrap function of fortran lapack routine dsytrf.
    static inline
    void dsytrf(char uplo, int n, ModuleBase::matrix& a, int lda, int *ipiv,
                double *work, int lwork , int info)
    {
        double * aux = LapackConnector::transpose_matrix(a, n, lda) ;
        dsytrf_(&uplo, &n, aux, &lda, ipiv, work, &lwork, &info);
        LapackConnector::transpose_matrix(aux, a, n, lda);
    }

    // wrap function of fortran lapack routine dsytri.
    static inline
    void dsytri(char uplo, int n, ModuleBase::matrix& a, int lda, int *iwork,double * work, int info)
    {
        double * aux = LapackConnector::transpose_matrix(a, n, lda) ;
        dsytri_(&uplo, &n, aux, &lda, iwork, work, &info);
        LapackConnector::transpose_matrix(aux, a, n, lda);
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
    void dgesvd(const char jobu,const char jobvt,const int m,const int n,
                ModuleBase::matrix &a,const int lda,double *s,ModuleBase::matrix &u,const int ldu,
                ModuleBase::matrix &vt,const int ldvt,double *work,const int lwork,int info)
    {
        //Transpose the matrix to the fortran-form

        double *aux = LapackConnector::transpose_matrix(a, n, lda);
        double *uux = LapackConnector::transpose_matrix(u, m, ldu);
        double *vtux = LapackConnector::transpose_matrix(vt, n, ldvt);

        dgesvd_(&jobu, &jobvt , &m, &n, aux, &lda, s, uux, &ldu, vtux, &ldvt, work, &lwork, &info);

        LapackConnector::transpose_matrix(aux, a, n, lda);
        LapackConnector::transpose_matrix(uux, u, m, ldu);
        LapackConnector::transpose_matrix(vtux, vt, n, ldvt);

        delete[] aux;
        delete[] uux;
        delete[] vtux;

    }

    static inline
    void zpotrf(char uplo,int n,ModuleBase::ComplexMatrix &a,const int lda,int* info)
    {
        std::complex<double> *aux = LapackConnector::transpose(a, n, lda);
        zpotrf_( &uplo, &n, aux, &lda, info);
        LapackConnector::transpose(aux, a, n, lda);
        delete[] aux;
		return;
    }

    static inline
    void zpotri(char uplo,int n,ModuleBase::ComplexMatrix &a,const int lda,int* info)
    {
        std::complex<double> *aux = LapackConnector::transpose(a, n, lda);
        zpotri_( &uplo, &n, aux, &lda, info);
        LapackConnector::transpose(aux, a, n, lda);
        delete[] aux;
		return;
    }

    static inline
    void spotrf(char uplo,int n,ModuleBase::matrix &a, int lda,int *info)
	{
		double *aux = LapackConnector::transpose_matrix(a, n, lda);
		for(int i=0; i<4; i++)std::cout << "\n aux=" << aux[i]; 
		spotrf_( &uplo, &n, aux, &lda, info);
		for(int i=0; i<4; i++)std::cout << "\n aux=" << aux[i]; 
		LapackConnector::transpose_matrix(aux, a, n, lda);
		delete[] aux;
		return;
	}

	static inline
	void spotri(char uplo,int n,ModuleBase::matrix &a, int lda, int *info)
	{
		double *aux = LapackConnector::transpose_matrix(a, n, lda);
		for(int i=0; i<4; i++)std::cout << "\n aux=" << aux[i]; 
		spotri_( &uplo, &n, aux, &lda, info);
		for(int i=0; i<4; i++)std::cout << "\n aux=" << aux[i]; 
		LapackConnector::transpose_matrix(aux, a, n, lda);
		delete[] aux;
		return;
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
    static inline
    void pdgetrf(int m, int n, ModuleBase::matrix &a,int ia, int ja, int *desca, int *ipiv, int *info)
    {
	double *aux = LapackConnector::transpose_matrix(a, n, m);
	pdgetrf_( &m, &n, aux, &ia, &ja, desca, ipiv, info);
	LapackConnector::transpose_matrix(aux, a, n, m);
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

	static inline
	double ddot_real(const int &dim, const std::complex<double>* psi_L,
		const std::complex<double>* psi_R)
	{
    	double *pL,*pR;
    	pL=(double *)psi_L;
    	pR=(double *)psi_R;
    	return dot(2*dim,pL,1,pR,1);
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
