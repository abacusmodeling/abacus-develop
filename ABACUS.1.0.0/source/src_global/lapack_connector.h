// =============================================================================
//                          C++ Header File
// Project:         LapackConnector
// File:            LapackConnector.hpp
// Author:          sltk
// Comment:         LapackConnector provide the connector to the fortran Lapack routine.
// Warning:
// Start time:      2007-03-08
// Last modified:   2008-08-12 ywcui : add zhegvx
// 					2008-08-13 mohan : find bug,test.
// 					2008-09-03 mohan : Add zgesv
// 					2009-03-08 mohan : add ilaenv
//					2010-01-22 spshu : add dgesvd
// =============================================================================

#ifndef LAPACKCONNECTOR_HPP
#define LAPACKCONNECTOR_HPP

#include <new>
#include <stdexcept>
#include <iostream>
#include "complexmatrix.h"
#include "matrix.h"


extern "C"
{
    int ilaenv_(int* ispec,const char* name,const char* opts,
    const int* n1,const int* n2,const int* n3,const int* n4);
    void zgesv_(const int* n,const int* nrhs,const complex<double> *A,const int *lda,
                const int *ipiv,const complex<double> *B,const int* ldb,int* info);
    void zhegv_(const int* itype,const char* jobz,const char* uplo,const int* n,
                complex<double>* a,const int* lda,complex<double>* b,const int* ldb,
                double* w,complex<double>* work,int* lwork,double* rwork,int* info);
	// mohan add dsygv 2010-03-20, to compute the eigenfunctions and eigenvalues of a real symmetry matrix.
	void dsygv_(const int* itype, const char* jobz,const char* uplo, const int* n,
				double* a,const int* lda,double* b,const int* ldb,
	 			double* w,double* work,int* lwork,int* info);
	// mohan add sspgvx 2010-03-21, to compute the selected eigenvalues of a generalized symmetric/Hermitian 
	// definite generalized problems. 
	void sspgvx_(const int* itype, const char* jobz,const char* range,const char* uplo,
				const int* n,double *ap,double *bp,const double* vl,
				const int* vu, const int* il, const int* iu,const double* abstol,
				const int* m,double* w,double *z,const int *ldz,
				double *work,int* iwork,int* ifail,int* info);
    void zheev_(const char* jobz,const char* uplo,const int* n,complex<double> *a,
                const int* lda,double* w,complex<double >* work,const int* lwork,
                double* rwork,int* info);
    void dsytrf_(const char* uplo, const int* n, double * a, const int* lda,
                 int *ipiv,double *work, int* lwork ,int *info);
    void dsytri_(const char* uplo,const int* n,double *a, const int *lda,
                 int *ipiv, double * work,int *info);
    void dgesvd_(const char* jobu,const char* jobvt,const int *m,const int *n,double *a,
                 const int *lda,double *s,double *u,const int *ldu,double *vt,
                 const int *ldvt,double *work,const int *lwork,int *info);
    void zhegvx_(const int* itype,const char* jobz,const char* range,const char* uplo,
                 const int* n,complex<double> *a,const int* lda,complex<double> *b,
                 const int* ldb,const double* vl,const double* vu,const int* il,
                 const int* iu,const double* abstol,const int* m,double* w,
                 complex<double> *z,const int *ldz,complex<double> *work,const int* lwork,
                 double* rwork,int* iwork,int* ifail,int* info);
    // mohan add dsptrf and dsptri 2010-02-03, to compute inverse real symmetry indefinit matrix.
    void spotrf_(const char* uplo,const int* n, double* A, int* lda, int *info);
    void spotri_(const char* uplo,const int* n, double* A, int* lda, int *info);
    // Peize Lin add dsptrf and dsptri 2016-06-21, to compute inverse real symmetry indefinit matrix.
    void dpotrf_(const char* uplo,const int* n, double* A, const int* lda, int *info);
    void dpotri_(const char* uplo,const int* n, double* A, const int* lda, int *info);
    // mohan add zpotrf and zpotri 2010-02-03, to compute inverse complex hermit matrix.
    void zpotrf_(const char* uplo,const int* n,const complex<double> *A,const int* lda,const int* info);
    void zpotri_(const char* uplo,const int* n,const complex<double> *A,const int* lda,const int* info);
	// Peize Lin add ?axpy 2016-08-04, to compute y=ax+y
	void saxpy_(int *N, float *alpha, float *X, int *incX, float *Y, int *incY);
	void daxpy_(int *N, double *alpha, double *X, int *incX, double *Y, int *incY);
	void caxpy_(int *N, complex<float> *alpha, complex<float> *X, int *incX, complex<float> *Y, int *incY);
	void zaxpy_(int *N, complex<double> *alpha, complex<double> *X, int *incX, complex<double> *Y, int *incY);
	// Peize Lin add ?scal 2016-08-04, to compute x=ax
	void sscal_(int *N, float *alpha, float *X, int *incX);
	void dscal_(int *N, double *alpha, double *X, int *incX);
	void cscal_(int *N, complex<float> *alpha, complex<float> *X, int *incX);
	void zscal_(int *N, complex<double> *alpha, complex<double> *X, int *incX);
};

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
class LapackConnector
{
private:
    // Transpose the complex matrix to the fortran-form real-complex array.
    static inline
    complex<double>* transpose(const ComplexMatrix& a, const int n, const int lda)
    {
        complex<double>* aux = new complex<double>[lda*n];
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
    double * transpose_matrix(const matrix& a, const int n, const int lda)
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

    // Transpose the fortran-form real-complex array to the complex matrix.
    static inline
    void transpose(const complex<double>* aux, ComplexMatrix& a, const int n, const int lda)
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
    void transpose_matrix(const double *aux, matrix &a, const int n, const int lda)
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
			default: throw invalid_argument("Side must be 'L' or 'R'");
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
			default: throw invalid_argument("Uplo must be 'U' or 'L'");
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
    void zgesv (const int n,const int nrhs,ComplexMatrix &A,const int lda,
                const int *ipiv,complex<double> *B,const int ldb,int* info)
    {
        complex<double> *aux = LapackConnector::transpose(A, n, lda);
        zgesv_(&n, &nrhs, aux, &lda, ipiv, B, &ldb, info);
        LapackConnector::transpose(aux, A, n, lda);
    }

    // wrap function of fortran lapack routine zhegv.
    static inline
    void zhegv(	const int itype,const char jobz,const char uplo,const int n,ComplexMatrix& a,
                const int lda,ComplexMatrix& b,const int ldb,double* w,complex<double>* work,
                int lwork,double* rwork,int info)
    {	// Transpose the complex matrix to the fortran-form real-complex array.
        complex<double>* aux = LapackConnector::transpose(a, n, lda);
        complex<double>* bux = LapackConnector::transpose(b, n, ldb);

        // call the fortran routine
        zhegv_(&itype, &jobz, &uplo, &n, aux, &lda, bux, &ldb, w, work, &lwork, rwork, &info);
        // Transpose the fortran-form real-complex array to the complex matrix.
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
				const int n,const matrix &ap,const matrix &bp,const double vl,
				const int vu, const int il, const int iu,const double abstol,
				const int m,double* w,matrix &z,const int ldz,
				double *work,int* iwork,int* ifail,int* info)
	{
	    // Transpose the complex matrix to the fortran-form real*16 array.
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
    void dsygv(	const int itype,const char jobz,const char uplo,const int n,matrix& a,
                const int lda,matrix& b,const int ldb,double* w,double* work,
                int lwork,int *info)
    {	// Transpose the complex matrix to the fortran-form real-complex array.
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


    // wrap function of fortran lapack routine zheev.
    static inline
    void zheev( const char jobz,
                const char uplo,
                const int n,
                ComplexMatrix& a,
                const int lda,
                double* w,
                complex< double >* work,
                const int lwork,
                double* rwork,
                int *info	)
    {	// Transpose the complex matrix to the fortran-form real-complex array.
        complex<double> *aux = LapackConnector::transpose(a, n, lda);
        // call the fortran routine
        zheev_(&jobz, &uplo, &n, aux, &lda, w, work, &lwork, rwork, info);
        // Transpose the fortran-form real-complex array to the complex matrix.
        LapackConnector::transpose(aux, a, n, lda);
        // free the memory.
        delete[] aux;
    }
	
    // wrap function of fortran lapack routine dsytrf.
    static inline
    void dsytrf(char uplo, int n, matrix& a, int lda, int *ipiv,
                double *work, int lwork , int info)
    {
        double * aux = LapackConnector::transpose_matrix(a, n, lda) ;
        dsytrf_(&uplo, &n, aux, &lda, ipiv, work, &lwork, &info);
        LapackConnector::transpose_matrix(aux, a, n, lda);
    }

    // wrap function of fortran lapack routine dsytri.
    static inline
    void dsytri(char uplo, int n, matrix& a, int lda, int *iwork,double * work, int info)
    {
        double * aux = LapackConnector::transpose_matrix(a, n, lda) ;
        dsytri_(&uplo, &n, aux, &lda, iwork, work, &info);
        LapackConnector::transpose_matrix(aux, a, n, lda);
    }


    // wrap function of fortran lapack routine zheev.
    static inline
    void zhegvx( const int itype, const char jobz, const char range, const char uplo,
                 const int n, const ComplexMatrix& a, const int lda, const ComplexMatrix& b,
                 const int ldb, const double vl, const double vu, const int il, const int iu,
                 const double abstol, const int m, double* w, ComplexMatrix& z, const int ldz,
                 complex<double>* work, const int lwork, double* rwork, int* iwork,
                 int* ifail, int info)
    {
        // Transpose the complex matrix to the fortran-form real-complex array.
        complex<double>* aux = LapackConnector::transpose(a, n, lda);
        complex<double>* bux = LapackConnector::transpose(b, n, ldb);
        complex<double>* zux = new complex<double>[n*iu];// mohan modify 2009-08-02
        //for(int i=0; i<n*iu; i++) zux[i] = complex<double>(0.0,0.0);

        // call the fortran routine
        zhegvx_(&itype, &jobz, &range, &uplo, &n, aux, &lda, bux, &ldb, &vl,
                &vu, &il,&iu, &abstol, &m, w, zux, &ldz, work, &lwork, rwork, iwork, ifail, &info);

        // Transpose the fortran-form real-complex array to the complex matrix
        LapackConnector::transpose(zux, z, iu, n);	 // mohan modify 2009-08-02

        // free the memory.
        delete[] aux;
        delete[] bux;
        delete[] zux;

    }
    static inline
    void dgesvd(const char jobu,const char jobvt,const int m,const int n,
                matrix &a,const int lda,double *s,matrix &u,const int ldu,
                matrix &vt,const int ldvt,double *work,const int lwork,int info)
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
    void zpotrf(char uplo,int n,ComplexMatrix &a,const int lda,int* info)
    {
        complex<double> *aux = LapackConnector::transpose(a, n, lda);
        zpotrf_( &uplo, &n, aux, &lda, info);
        LapackConnector::transpose(aux, a, n, lda);
        delete[] aux;
		return;
    }

    static inline
    void zpotri(char uplo,int n,ComplexMatrix &a,const int lda,int* info)
    {
        complex<double> *aux = LapackConnector::transpose(a, n, lda);
        zpotri_( &uplo, &n, aux, &lda, info);
        LapackConnector::transpose(aux, a, n, lda);
        delete[] aux;
		return;
    }

    static inline
    void spotrf(char uplo,int n,matrix &a, int lda,int *info)
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
	void spotri(char uplo,int n,matrix &a, int lda, int *info)
	{
		double *aux = LapackConnector::transpose_matrix(a, n, lda);
		for(int i=0; i<4; i++)std::cout << "\n aux=" << aux[i]; 
		spotri_( &uplo, &n, aux, &lda, info);
		for(int i=0; i<4; i++)std::cout << "\n aux=" << aux[i]; 
		LapackConnector::transpose_matrix(aux, a, n, lda);
		delete[] aux;
		return;
	}

	// Peize Lin add 2016-07-09
	static inline
	void dpotrf( char uplo, const int n, matrix &a, const int lda, int *info )
	{
		const char uplo_changed = change_uplo(uplo);
		dpotrf_( &uplo_changed, &n, a.c, &lda, info );
	}	
	
	// Peize Lin add 2016-07-09
	static inline
	void dpotri( char uplo, const int n, matrix &a, const int lda, int *info )
	{
		const char uplo_changed = change_uplo(uplo);
		dpotri_( &uplo_changed, &n, a.c, &lda, info);		
	}	
	
	// Peize Lin add 2016-08-04
	static inline 
	void axpy( int N,  float alpha,  float *X,  int incX, float *Y,  int incY)
	{
		saxpy_(&N, &alpha, X, &incX, Y, &incY);
	}	
	static inline 
	void axpy( int N,  double alpha,  double *X,  int incX, double *Y,  int incY)
	{
		daxpy_(&N, &alpha, X, &incX, Y, &incY);
	}	
	static inline 
	void axpy( int N,  complex<float> alpha,  complex<float> *X,  int incX, complex<float> *Y,  int incY)
	{
		caxpy_(&N, &alpha, X, &incX, Y, &incY);
	}	
	static inline 
	void axpy( int N,  complex<double> alpha,  complex<double> *X,  int incX, complex<double> *Y,  int incY)
	{
		zaxpy_(&N, &alpha, X, &incX, Y, &incY);
	}	
	
	// Peize Lin add 2016-08-04
	static inline 
	void scal( int N,  float alpha, float *X,  int incX)
	{
		sscal_(&N, &alpha, X, &incX);
	}	
	static inline 
	void scal( int N,  double alpha, double *X,  int incX)
	{
		dscal_(&N, &alpha, X, &incX);
	}	
	static inline 
	void scal( int N,  complex<float> alpha, complex<float> *X,  int incX)
	{
		cscal_(&N, &alpha, X, &incX);
	}	
	static inline 
	void scal( int N,  complex<double> alpha, complex<double> *X,  int incX)
	{
		zscal_(&N, &alpha, X, &incX);
	}	

};
#endif  // LAPACKCONNECTOR_HPP
