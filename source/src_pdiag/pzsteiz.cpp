#include "pzsteiz.h"

void pzsteiz(int n,double *a,double *b,double *w,std::complex <double> *Z,int m)
/*
 * PSEPS routine (version 2.0) --
 * Computer Network Information Center, CAS. 
 * July 15, 2006
 * 
 * Purpose
 *========
 * pzsteiz computes the eigenvector of a symmetric tridiagonal matrix .
 * Basically a parallel version of LAPACK's ZSTEIN
 *
 * Arguments
 * ========= 
 * n (global input) int  
 *   The order of the tridiagonal matrix T.  n >= 0.  
 * a (global input) double array, dimension (n)
 *   The n diagonal elements of the tridiagonal matrix T.
 * b (global input) double precision array, dimension (n-1)
 *   The (n-1) off-diagonal elements of the tridiagonal matrix T.
 * w (local input) double array, dim (m)
 *   the first m elements of w contain the eigenvalues for which eigenvectors are to be computed.
 * Z  (local output) double precision std::complex array, dimension (m*n)
 *    Z contains the computed eigenvectors associated with the  specified eigenvalues. 
 */
{
    TITLE("Parallel","pzsteiz");
    int i,j,info,ldz=n;

    double work[5*n];
    int iblock[n],iwork[n],ifail[n],isplit[n];
    info=0;
    isplit[0]=n;
    for (i=0; i<n; i++)
	{
        iblock[i]=1;
	}
    zstein_(&n,a,b,&m,w,iblock,isplit,Z,&ldz,work,iwork,ifail,&info);
}
