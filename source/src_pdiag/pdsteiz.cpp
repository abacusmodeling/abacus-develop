#include "pdsteiz.h"

void pdsteiz(int n,double *a,double *b,double *w,double *Z,int m)
/*
 * PSEPS routine (version 1.0) --
 * Computer Network Information Center, CAS.
 * December 15, 2004
 *
 * Purpose
 *========
 * psteiz computes the eigenvector of a symmetric tridiagonal matrix .
 * Basically a parallel version of LAPACK's DSTEIN
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
 * Z (local output) double precision array, dimension (m*n)
 *   Z contains the computed eigenvectors associated with the  specified eigenvalues.
 */
{
	TITLE("Parallel_Diago","pdsteiz");
    int i,j,info,ldz=n;

    double work[5*n];
    int iblock[n],iwork[n],ifail[m],isplit[n];
    info=0;
    isplit[0]=n;
    for (i=0; i<m; i++)
	{
        iblock[i]=1;
	}
    dstein_(&n,a,b,&m,w,iblock,isplit,Z,&ldz,work,iwork,ifail,&info);
	return;
}
