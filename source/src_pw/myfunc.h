#ifndef MYFUNC_H
#define MYFUNC_H
using namespace std;

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/vector3.h"
#include "../module_base/complexmatrix.h"
//#include "../src_algorithms/mymath.h"     //  only wgauss(),wlgauss()
#include <complex>

void gcxc(double rho, double grho, double &sx, double &sc,
          double &v1x, double &v2x, double &v1c, double &v2c);
void becke88(double rho, double grho, double &sx, double &v1x, double &v2x);
void perdew86(double rho, double grho, double &sc, double &v1c, double &v2c);
void glyp(double rho, double grho, double &sc, double &v1c, double &v2c);
void hcth(double rho, double grho, double &sx, double &v1x, double &v2x);
void pwcorr(double r, double c[], double &g, double &dg);
void optx(double rho, double grho, double &sx, double &v1x, double &v2x);

// Consider spin :
void gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
              double &sx, double &v1xup, double &v1xdw, double &v2xup,
              double &v2xdw);
void gcc_spin(double rho, double &zeta, double grho, double &sc,
              double &v1cup, double &v1cdw, double &v2c);
void becke88_spin(double rho, double grho, double &sx, double &v1x,
                  double &v2x);
void perdew86_spin(double rho, double zeta, double grho, double &sc,
                   double &v1cup, double &v1cdw, double &v2c);
void ggac_spin(double rho, double zeta, double grho, double &sc,
               double &v1cup, double &v1cdw, double &v2c);
void pbec_spin(double rho, double zeta, double grho, const int &flag, double &sc,
               double &v1cup, double &v1cdw, double &v2c);

// GGA : gradcorr (myfunc3.cpp)
void gradcorr(double &etxc, double &vtxc, ModuleBase::matrix &v);

void gradient(int ncx, int ncy, int ncz, int ncxyz, double *a, int ngm,
              ModuleBase::Vector3 < double> *g, int *nl, double lat0,
              ModuleBase::Vector3 < double> *ga);
//void grad_dot (int ncx, int ncy, int ncz,int ncxyz, ModuleBase::matrix a, int ngm,
//	       ModuleBase::Vector3 < double> *g, int *nl, double lat0,double *da);
void grad_dot(int ncx, int ncy, int ncz, int ncxyz, ModuleBase::Vector3 < double> *a, int ngmc,
              ModuleBase::Vector3 < double> *g, int *ig2fftc, double lat0, double *da);


// myfunc5.cpp (from LPACK)
// compute y = alpha * y where alpha is a scalar and y is an n-std::vector
void dscal(const int n, const double &alpha, double *y, const int incy);
// a(i,:) = alpha * a(i,:) where a is a matrix
void dscal(const double &alpha, const ModuleBase::matrix &a, const int i);
// compute y = alpha * x + y where alpha is a scalar and x and y are n-vectors
void daxpy(const int n, const double &alpha, const double *x, const int incx, double *y, const int incy);
void zaxpy(int n, double alpha, std::complex < double> *x, int incx,std::complex<double> *y, int incy);
void zaxpy(int n, std::complex < double> alpha, std::complex < double> *x, int incx,std::complex<double> *y, int incy);
// y(i,:) = alpha * x + y(i,:) where y is a matrix
void zaxpy(double alpha,std::complex<double> *x,ModuleBase::ComplexMatrix &y,int i);
// y(i,:) = alpha * x + y(i,:)
void zaxpy(double alpha,const ModuleBase::ComplexMatrix &x,int i,std::complex<double> *y);
void zaxpy(std::complex<double> alpha,const ModuleBase::ComplexMatrix &x,int i,std::complex<double> *y);
// y(j,:) = alpha * x(i,:) + y(j,:)
void zaxpy(std::complex<double> alpha,const ModuleBase::ComplexMatrix &x,int i,ModuleBase::ComplexMatrix &y,int j);
// copy x to y, where x and y is n-vectors



//*********************************************************************************
//                             dcopy
//*********************************************************************************
void dcopy(int n,std::complex<double> *x,int incx,std::complex<double> *y,int incy);
void dcopy(int n,double *x,int incx,double *y,int incy);
void dcopy(int n,int *x,int incx,int *y,int incy);

// copy ith row of a to y where a is a matrix and y is a std::vector
void dcopy(const ModuleBase::ComplexMatrix &a,int i,std::complex<double> *y);
void dcopy(const ModuleBase::matrix &a,int i,double *y);
void dcopy(const ModuleBase::matrix &a,int i,int *y);

// copy x to ith row of b where b is a matrix and x is a std::vector
void dcopy(std::complex<double> *x,ModuleBase::ComplexMatrix &b,int i);
void dcopy(double *x,ModuleBase::matrix &b,int i);

// b(j,:) = a(i,:)
void dcopy(const ModuleBase::ComplexMatrix &a,int i,ModuleBase::ComplexMatrix &b,int j);
void dcopy(const ModuleBase::matrix &a,       int i,ModuleBase::matrix &b,       int j);
void dcopy(int n,const ModuleBase::ComplexMatrix &a,int inca,ModuleBase::ComplexMatrix &b,int i,int incb);



//*********************************************************************************
//                             ddot
//*********************************************************************************
double ddot(int n,double *x,int incx,double *y,int incy);
std::complex<double> ddot(int n,std::complex<double> *x,int incx,std::complex<double> *y, int incy);
std::complex<double> ddot(const ModuleBase::ComplexMatrix &a,int i,std::complex<double> *x);
std::complex<double> ddot(std::complex<double> *x,const ModuleBase::ComplexMatrix &y,int i);
std::complex<double> ddot(const ModuleBase::ComplexMatrix &x,int i,const ModuleBase::ComplexMatrix &y,int j);
double ddot(const ModuleBase::matrix &a,int i,double *y);
double ddot(const ModuleBase::matrix &x,int i,const ModuleBase::matrix &y,int j);


// perform one of the matrix-matrix operations
// C = alpha * op(A) * op(B) + beta * C
void dgemm(char tra, char trb, int m, int n, int k, double alpha,const ModuleBase::matrix A, int lda, const ModuleBase::matrix B, int ldb, double beta,ModuleBase::matrix &C, int ldc);
void zgemm(char tra,char trb,int m,int n,int k,std::complex<double> alpha,
           const ModuleBase::ComplexMatrix &A,
           int lda,
           const ModuleBase::ComplexMatrix &B,
           int ldb,
           std::complex < double> beta,
           ModuleBase::ComplexMatrix &c,
           int ldc);
void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           std::complex < double> alpha,
           const std::complex<double> *A,
           int lda,
           const ModuleBase::ComplexMatrix &B,
           int ldb,
           std::complex < double> beta,
           ModuleBase::ComplexMatrix &c,
           int ldc);

void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           std::complex < double> alpha,
           const std::complex<double> *A,
           int lda,
           const ModuleBase::ComplexMatrix &B,
           int ldb,
           std::complex < double> beta,
           std::complex <double> *c,
           int ldc);

void zgemv(char , int , int , std::complex < double> alpha ,
           ModuleBase::ComplexMatrix overlap, int , std::complex < double> swfcatom , int npwx,
           std::complex < double>  , ModuleBase::ComplexMatrix work, int);//called in orthoatwfc()

int ILAENV(int , char *name, char *opts,const int n1,const int n2,const int n3,const int n4);

void ZHPEV(int , std::complex < double> *hp, double *e, ModuleBase::ComplexMatrix &v,
           int ldh, int n, std::complex < double> *aux, int naux);

// compute the Euclidean length (12 norm) of std::vector x, with scaling of
// input to avoid destructive underflow and overflow
double dnrm2(const int n, const double *x, const int incx) ;
double dnrm2(const ModuleBase::matrix &a, int i) ;	// i-row of matrix a

// additional
//---------------------------------

double min3(double x1, double x2, double x3);

// unsolved :
//---------------------------------------------------------------
// for US-PP only
void addusdens();

// LDA+U only (not used here)
void new_ns();


// symmetry (not used here)
void psymrho(ModuleBase::matrix rho, int nvx, int ncy, int ncz, int nsym, double s,
             int ftau);
void symrho(ModuleBase::matrix rho, int ncx, int ncy, int ncz, int nsym, double s,
            int ftau);

// Parallel code only
void poolrecover(ModuleBase::matrix et, int nbnd, int nkstot, int nks);
void poolreduce(int, double eband);
void poolscatter(int nbnd, int nkstot, ModuleBase::matrix , int nks, ModuleBase::matrix wg);
void mp_bcast(double ef, int root_image, int intra_image_comm);
void mp_bcast(double *e, int root_image, int intra_image_comm);
void mp_bcast(ModuleBase::ComplexMatrix, int root_image, int intra_image_comm);

// not used
void vhpsi(int lda, int n, int m, ModuleBase::ComplexMatrix psi, ModuleBase::ComplexMatrix hpsi);

// Later
void c_phase();
void ns_adj();
void vpack(int NCXYZ, int ncxyz, int nspin, double *vnew, double *vr, int);

// here is not definition
//-------------------------------
//void REWIND( int iunigk );

void inquire(int unit, int opened);
void divide(int nqxq, int startq, int lastq);
void reduce(int, double *dr2);
void ireduce(int , int ngkp);
void davcio(ModuleBase::ComplexMatrix evc, int nwordwfc, int iunwfc, int ik, int);

void ZHEGVX(
    int itype,
    char jobz,
    char range ,
    char uplo ,
    const int n,
    const ModuleBase::ComplexMatrix &a,
    const int lda,
    const ModuleBase::ComplexMatrix &b,
    const int ldb,
    double vl,
    double vu,
    int il ,
    int iu,
    double abstol,
    int &m,
    double *w,
    ModuleBase::ComplexMatrix &z,
    const int ldz,
    double *work,
    int lwork,
    double *rwork,
    int *iwork,
    int *ifail,
    int &info );


#endif // NYFUNC
void DGER(int na, int nb, double , ModuleBase::ComplexMatrix a, int lda,
          ModuleBase::ComplexMatrix b, int ldb, ModuleBase::ComplexMatrix c, int ldc);

void add_efield(ModuleBase::matrix rho, double *v, double etotefield);	// no used



