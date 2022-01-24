/* myfunc.cpp */
// from LPACK
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

using namespace std;
#include "myfunc.h"
#include "../module_base/blas_connector.h"
#include "global.h"

// dscal compute y = alpha * y, where alpha is a scalar and
void dscal(const int n,const double &alpha,double *y,const int incy)
{
    if (incy <= 0 || n < 1)
    {
        std::cout << "\n error in dscal,incy <= 0 or n < 1";
        return;
    }

    for (int i = 0; i < n; i += incy)
    {
        y[i] *= alpha;
    }

    return;
}

// a(i,:) = alpha * a(i,:) where a is a matrix
// i line
void dscal(const double &alpha,ModuleBase::matrix &a,const int i)
{
    int nc = a.nc;
    int nr = a.nr;

    if (nc <= 0 || nr <= 0 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dscal,nc <= 0 or nr <= 0 or i < 0 or i >= nr, ";
        return;
    }

    for (int j = 0;j < nc;j++)
    {
        a(i,j) = a(i,j) * alpha;
    }
}

// daxpy compute y := alpha * x + y where alpha is a scalar and x and y
void daxpy(const int n, const double &alpha, const double *x, const int incx, double *y, const int incy)
{
    if (n < 1 || incy <= 0 || incx <= 0)
    {
        std::cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }
    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] += alpha * x[ix];
    }
    return;
}

void zaxpy(int n, double alpha, std::complex < double> *x, int incx, std::complex < double> *y, int incy)
{

    if (n < 1 || incy <= 0 || incx <= 0)
    {
        std::cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)  // y := alpha * x + y
    {
        y[iy] += alpha * x[ix];
    }
}


//-------------------------------------------------------------
void zaxpy(int n, std::complex < double> alpha,  std::complex < double> *x,
           int incx, std::complex < double> *y, int incy)
{
    // zaxpy compute y := alpha * x + y where alpha is a scalar and
    // x and y are n-vectors.
    // DOUBLE COMPLEX ALPHA
    // DOUBLE COMPLEX X(*), Y(*)
    // INTEGER N, INCX, INCY

    // ARGUMENTS

    // N (input)
    //		On entry, N specifies the number  of  elements  in
    //		the  std::vector.   N must be at least one for the sub-
    //		routine to have any visible effect.  Unchanged  on
    //		exit.

    // ALPHA (input)
    //		On  entry,  ALPHA  specifies  the  scalar   alpha.
    //		Unchanged on exit.

    // X (input)
    //		array of DIMENSION at least ( 1 + ( n -  1  )*abs(
    //		INCX  )  ).  Before entry, the incremented array X
    //		must contain the std::vector x.  Unchanged on exit.

    // INCX (input)
    //		On entry, INCX specifies  the  increment  for  the
    //		elements of X. Unchanged on exit.

    // Y (input/output)
    //		array of DIMENSION at least ( 1 + ( n -  1  )*abs(
    //		INCY  ) ).  On entry, the incremented array Y must
    //		contain the std::vector y. On exit, Y is overwritten by
    //		the updated std::vector y.

    // INCY (input)
    //		On entry, INCY specifies  the  increment  for  the
    //		elements of Y. Unchanged on exit.

    if (n < 1 || incy <= 0 || incx <= 0)
    {
        std::cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)  // y := alpha * x + y
    {
        y[iy] += alpha * x[ix];
    }

} // end zaxpy

// y(i,:) = alpha * x + y(i,:) where y is a matrix
void zaxpy(double alpha,
           std::complex < double> *x,
           ModuleBase::ComplexMatrix &y,
           int i)
{
    int nr, nc;
    nr = y.nr;
    nc = y.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in daxpy, nr < 1 or nc < 1 or i out of range, ";
        return;
    }

    for (int j = 0;j < nc;j ++) // y := alpha * x + y
    {
        y(i, j) += alpha * x[j];
    }

//    std::cout << "\n End daxpy() " << std::endl;
}

// y = alpha * x(i,:) + y
void zaxpy(double alpha,
           const ModuleBase::ComplexMatrix &x,
           int i,
           std::complex < double> *y)
{
    int nr, nc;
    nr = x.nr;
    nc = x.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in daxpy, nr < 1 or nc < 1 or i out of range, ";
        return;
    }

    for (int j = 0;j < nc;j ++) // y := alpha * x + y
    {
        y[j] += alpha * x(i, j);
    }
}

void zaxpy(std::complex < double> alpha,
           const ModuleBase::ComplexMatrix &x,
           int i,
           std::complex < double> *y)
{
    int nr, nc;
    nr = x.nr;
    nc = x.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in daxpy, nr < 1 or nc < 1 or i out of range, ";
        return;
    }

    for (int j = 0;j < nc;j ++) // y := alpha * x + y
    {
        y[j] += alpha * x(i, j);
    }
}

// y(j,:) = alpha * x(i,:) + y(j,:)
void zaxpy(std::complex < double> alpha,
           const ModuleBase::ComplexMatrix &x,
           int i,
           ModuleBase::ComplexMatrix &y,
           int j)
{
    int nr, nc;
    nr = y.nr;
    nc = y.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr || j < 0 || j >= nr)
    {
        std::cerr << "\n error in daxpy, nr < 1 or nc < 1 or i or j out of range, ";
        return;
    }

    for (int k = 0;k < nc;k ++) // y := alpha * x + y
    {
        y(j, k) += alpha * x(i, k);
    }
}

//-----------------------------------------------------------------
void dcopy(int n, double *x, int incx, double *y, int incy)
{
    // dcopy Copy x to y where x and y are n-vectors.
    if (n < 1 || incx <= 0 || incy <= 0)
    {
        std::cerr << "\n error in dcopy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] = x[ix];
    }
} // end dcopy

//------------------------------------------------------------------
void dcopy(int n, std::complex < double> *x, int incx, std::complex < double> *y, int incy)
{
    // zcopy Copy x to y where x and y are n-vectors.
    if (n < 1 || incx <= 0 || incy <= 0)
    {
        std::cerr << "\n error in dcopy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] = x[ix];
    }

} // end dcopy

void dcopy(int n, int *x, int incx, int *y, int incy)
{
    if (n < 1 || incx <= 0 || incy <= 0)
    {
        std::cerr << "\n error in dcopy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] = x[ix];
    }
}

/* Copy a(i,:) to y where x is matrix, and y are n-vectors. */
void dcopy(const ModuleBase::matrix &a,
           int i,
           double *y)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy((matrix a, int i, double *y)), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = a(i, j);
    }

} // end dcopy

//-------------------------------------
void dcopy(const ModuleBase::matrix &a,
           int i,
           int *y)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(matrix a, int i, int *y),"
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = (int) a(i, j);
    }

} // end dcopy

//-----------------------------------------------------------------
void dcopy(const ModuleBase::ComplexMatrix &a,
           int i,
           std::complex < double> *y)
{
    // dcopy Copy a(i,:) to y where a is std::complex matrix, and y are n-vectors.
    const int nr = a.nr;
    const int nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(ModuleBase::ComplexMatrix a, int i, std::complex < double> *),"
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = a(i, j);
    }
} // end dcopy

// ------------------------------------
void dcopy(double *x, ModuleBase::matrix &b, int i)
{
    // copy x to ith row of b where b is a matrix and x is a std::vector
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(double *x, matrix &b, int i), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        b(i, j) = x[j];
    }
}

void dcopy(std::complex < double> *x, ModuleBase::ComplexMatrix &b, int i)
{
    // copy x to ith row of b where b is a std::complex matrix and x is a std::vector
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(double> *x, ModuleBase::ComplexMatrix &b, int i), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        b(i, j) = x[j];
    }
}

// b(j,:) = a(i,:)
void dcopy(const ModuleBase::matrix &a,
           int i,
           ModuleBase::matrix &b,
           int j)
{
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(matrix a, int i,matrix &b, int ), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int k = 0;k < nc; k++)
    {
        b(j, k) = a(i, k);
    }
}

void dcopy(const ModuleBase::ComplexMatrix &a,
           int i,
           ModuleBase::ComplexMatrix &b,
           int j)
{
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(ModuleBase::ComplexMatrix a, int i,ModuleBase::ComplexMatrix &b, int j), "
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int k = 0;k < nc; k++)
    {
        b(j, k) = a(i, k);
    }
}

//---------------------------------------------------------
double ddot(int n,
            double *x,
            int incx,
            double *y,
            int incy)
{
    // ddot compute the dot product of x and y where x  and  y  are
    // n-vectors.
    double prod;

    if (n < 1 || incx <= 0 || incy <= 0)
    {
        std::cerr << "\n error in ddot, n < 1 or incx <= 0 or incy <= 0, ";
        return 0;
    }

    prod = 0.0;

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        prod += x[ix] * y[iy];
    }

    return prod;
} // end ddot

//-----------------------------------------------------------------------
std::complex < double> ddot(int n,
                       std::complex < double> *x,
                       int incx,
                       std::complex < double> *y,
                       int incy)
{
    // zdotc compute the dot product of conjg(x) and y where x  and
    // y are n-vectors.
    std::complex < double> prod;

    if (n < 1 || incx <= 0 || incy <= 0)
    {
        std::cerr << "\n error in ddot, n < 1 or incx <= 0 or incy <= 0, ";
        return 0;
    }

    prod = std::complex < double>(0.0, 0.0);

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        prod += conj(x[ix]) * y[iy];
    }

    return prod;
} // end zdotc

// -----------------------------------------------------------------------
std::complex < double> ddot(const ModuleBase::ComplexMatrix &a,
                       int i,
                       std::complex < double> *y)
{
    //compute the dot product of i_th row of matrix a and y where
    // y are a std::vector.
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr <= 1 || nc <= 1 || i < 0 || i >= nr)
    {
        std::cout << "\n error in ddot, nr or nc < 1 or i out of range ";
        return 0;
    }

    std::complex < double> z;

    z = std::complex < double> (0, 0);

    for (int k = 0; k < nc; k++)
    {
        z += conj(a(i, k)) * y[k];
    }

    return z;
}

//--------------------------------
double ddot(const ModuleBase::matrix &a,
            int i,
            double  *y)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr <= 1 || nc <= 1 || i < 0 || i >= nr)
    {
        std::cout << "\n error in ddot, nr or nc < 1 or i out of range ";
        return 0;
    }

    double z = 0;

    for (int k = 0; k < nc; k++)
    {
        z += a(i, k) * y[k];
    }

    return z;
}

double ddot(const ModuleBase::matrix &a,
            int i,
            const ModuleBase::matrix &b,
            int j)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr <= 1 || nc <= 1 || i < 0 || i >= nr ||
            j < 0 || j >= nr)
    {
        std::cout << "\n error in ddot, nr or nc < 1 or i or j out of range ";
        return 0;
    }

    double z = 0;

    for (int k = 0; k < nc; k++)
    {
        z += a(i, k) * b(j, k);
    }

    return z;
}

double dnrm2(const int n, const double *x, const int incx)
{
    // compute Euclidean length (12 norm) of std::vector x,
    if (n < 0 || incx <= 0)
    {
        std::cerr << "\n error in dnrm2, n < 0 or incx <= 0, ";
        return 0;
    }
    if (n == 0)
    {
        return 0;
    }

    double norm2=0.0;
    for (int ix=0; ix<n; ix+=incx)
    {
        norm2 += x[ix] * x[ix];
    }

    return sqrt(norm2);
}

// i-row of matrix a
double dnrm2(const ModuleBase::matrix &a,
             int i)
{
    int nc, nr;
    nc = a.nc;
    nr = a.nr;
    double norm2;
    norm2 = 0.0;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dnrm2, nr or nc < 1 or i out of range ";
        return 0;
    }

    for (int j = 0;j < nc; j++)
    {
        norm2 += a(i, j) * a(i, j);
    }

    return sqrt(norm2);
}

void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           std::complex < double> alpha,
           const std::complex<double> *a,
           int lda,
           const ModuleBase::ComplexMatrix &b,
           int ldb,
           std::complex < double> beta,
           std::complex<double> *c,
           int ldc)
{
	ModuleBase::TITLE("myfunc5","zgemm1");
    //int nra=1;
	int nca=k;
	int nrb=b.nr;
	int ncb=b.nc;
	int nrc=1;
	int ncc=n;

    //double alpha0[2] = {alpha.real(), alpha.imag()};
	//double beta0[2] = {beta.real(), beta.imag()};

    double *aux, *bux, *cux; //HLX: bug fixed on 12/20/2006
    aux = new double[2*nca];
    bux = new double[2*nrb*ncb];
    cux = new double[2*ncc];

    int i;
    int j;
    int ij;

    for (i = 0; i < nca; i++)
    {
        ij = 2*i;
        aux[ij] = a[i].real();
        aux[ij+1] = a[i].imag();
    }
    for (i = 0; i < ncb; i++)
    {
        for (j = 0; j < nrb; j++)
        {
            ij = 2 * (j + i * nrb);
            bux[ij]  = b(j, i).real();
            bux[ij+1] = b(j, i).imag();
        }
    }


    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * i;
            cux[ij]  = c[i].real();
            cux[ij+1] = c[i].imag();
        }
    }

	//ModuleBase::WARNING_QUIT("reset zgemm","reset zgemm");
//    zgemm_(&tra, &trb, &m, &n, &k, alpha0, aux, &lda, bux, &ldb, beta0, cux, &ldc);
	ModuleBase::WARNING_QUIT("myfunc_5::zgemm","please don't ever use it again.");

    for (i = 0; i < ncc; i++)
    {
        ij = 2 *  i * ldc;
        c[i] = cux[ij] + std::complex< double>(0, 1) * cux[ij+1];
    }

    delete [] aux;
    delete [] bux;
    delete [] cux;
}


void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           std::complex < double> alpha,
           const std::complex<double> *a,
           int lda,
           const ModuleBase::ComplexMatrix &b,
           int ldb,
           std::complex < double> beta,
           ModuleBase::ComplexMatrix &c,
           int ldc)
{
    ModuleBase::TITLE("myfunc5","zgemm2");
    //int nra = 1;
	int nca = k;
	int nrb = b.nr;
	int ncb = b.nc;
	int nrc = c.nr;
	int ncc = c.nc;

	//double alpha0[2] = {alpha.real(), alpha.imag()};
	//double beta0[2] = {beta.real(), beta.imag()};

    double *aux, *bux, *cux; //HLX: bug fixed on 12/20/2006
    aux = new double[2*nca];
    bux = new double[2*nrb*ncb];
    cux = new double[2*nrc*ncc];

    int i;
    int j;
    int ij;

    for (i = 0; i < nca; i++)
    {
        ij = 2*i;
        aux[ij] = a[i].real();
        aux[ij+1] = a[i].imag();
    }
    for (i = 0; i < ncb; i++)
    {
        for (j = 0; j < nrb; j++)
        {
            ij = 2 * (j + i * nrb);
            bux[ij]  = b(j, i).real();
            bux[ij+1] = b(j, i).imag();
        }
    }


    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * (j + i * nrc);
            cux[ij]  = c(j, i).real();
            cux[ij+1] = c(j, i).imag();
        }
    }

	//ModuleBase::WARNING_QUIT("reset zgemm","reset zgemm");
    //zgemm_(&tra, &trb, &m, &n, &k, alpha0, aux, &lda, bux, &ldb, beta0, cux, &ldc);
	ModuleBase::WARNING_QUIT("myfunc_5::zgemm","please don't ever use it again.");

    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2*(j+i*ldc);
            c(j,i) = cux[ij] + std::complex< double>(0,1) * cux[ij+1];
            //		std::cout<<std::setw(12)<<c(j,i);
        }
        //	std::cout<<std::endl;
    }

    delete [] aux;
    delete [] bux;
    delete [] cux;
}


// C = alpha * op(A) * op(B) + beta * C
void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           std::complex < double> alpha,
           const ModuleBase::ComplexMatrix &a,
           int lda,
           const ModuleBase::ComplexMatrix &b,
           int ldb,
           std::complex < double> beta,
           ModuleBase::ComplexMatrix &c,
           int ldc)
{
//	ModuleBase::TITLE("myfunc5","zgemm3");
    int nra, nca, nrb, ncb, nrc, ncc ;
    nra = a.nr,   nca = a.nc;
    nrb = b.nr,   ncb = b.nc;
    nrc = c.nr,   ncc = c.nc;

	//double alpha0[2] = {alpha.real(), alpha.imag()};
	//double beta0[2] = {beta.real(), beta.imag()};

    double *aux, *bux, *cux; //HLX: bug fixed on 12/20/2006
    aux = new double[2*nra*nca];
    bux = new double[2*nrb*ncb];
    cux = new double[2*nrc*ncc];

    int i;
    int j;
    int ij;

    for (i = 0; i < nca; i++)
    {
        for (j = 0; j < nra; j++)
        {
            ij = 2 * (j + i * nra);
            aux[ij] = a(j, i).real();
            aux[ij+1] = a(j, i).imag();
        }
    }

    for (i = 0; i < ncb; i++)
    {
        for (j = 0; j < nrb; j++)
        {
            ij = 2 * (j + i * nrb);
            bux[ij]  = b(j, i).real();
            bux[ij+1] = b(j, i).imag();
        }
    }


    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * (j + i * nrc);
            cux[ij]  = c(j, i).real();
            cux[ij+1] = c(j, i).imag();
        }
    }

	//ModuleBase::WARNING_QUIT("reset zgemm","reset zgemm");
    //zgemm_(&tra, &trb, &m, &n, &k, alpha0, aux, &lda, bux, &ldb, beta0, cux, &ldc);
	ModuleBase::WARNING_QUIT("myfunc_5::zgemm","please don't ever use it again.");

    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * (j + i * ldc);
            c(j, i) = cux[ij] + std::complex< double>(0, 1) * cux[ij+1];
        }
    }

    delete [] aux;

    delete [] bux;
    delete [] cux;
}//end zgemm


// C = alpha * op(A) * op(B) + beta * C
void dgemm(char tra, char trb, int m, int n, int k, double alpha,
           const ModuleBase::matrix  a, int lda, const ModuleBase::matrix b, int ldb, double beta,
           ModuleBase::matrix &c, int ldc)
{
//    std::cout << "\n === ZGEMM() ===" << std::endl;
    ModuleBase::matrix a1;

    if (tra == 'n' || tra == 'N')
        a1.create(a.nr, a.nc);
    else
        a1.create(a.nc, a.nr);

//    std::cout << "\n a1.nr = " << a1.nr
//	 << "   a1.nc = " << a1.nc << std::endl;

    c = beta * c;

    if (tra == 'n' || tra == 'N')
        a1 = alpha * a;
    else
        a1 = alpha * transpose(a);

    if (trb == 'n' || trb == 'N')
        c += a1 * b;
    else
        c += a1 * transpose(b);

//    a1.freemem();
//    std::cout << "\n end ZGEMM() " << std::endl;
    return;
}//end dgemm

