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

// compute y = alpha * x + y where alpha is a scalar and x and y are n-vectors
void daxpy(const int n, const double &alpha, const double *x, const int incx, double *y, const int incy);
// copy ith row of a to y where a is a matrix and y is a std::vector
void dcopy(const ModuleBase::ComplexMatrix &a,int i,std::complex<double> *y);
// copy x to ith row of b where b is a matrix and x is a std::vector
void dcopy(double *x,ModuleBase::matrix &b,int i);

// compute the Euclidean length (12 norm) of std::vector x, with scaling of
// input to avoid destructive underflow and overflow
double dnrm2(const int n, const double *x, const int incx) ;
#endif // NYFUNC



