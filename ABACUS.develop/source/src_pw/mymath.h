/* mymath.h file */

#ifndef MYMATH_H
#define MYMATH_H
#include "../src_pw/tools.h"

#if defined __FFTW2
#include "fftw.h"
#elif defined __FFTW3
#include "fftw3.h"
#else
#include <fftw3-mpi.h>
//#include "fftw3-mpi_mkl.h"
#endif

// in mymath1.cpp
#ifdef __FFTW3
void fftw_zeros(fftw_complex *data,int n);
#endif


// in mymath1.cpp
double rndm();
void simpson(const int mesh,const double *func,const double *rab, double &asum);
//void simpson_cp90( int mesh, double *func, double *rab, double intg );
//void simpson_fpmd(int n, double *func, double dx, double s);

double my_erf(double x);
double my_erfc(double x);
double gauss_freq(double x);

// in mymath3.cpp
void heapsort(int n, double *r, int *ind);
void heapAjust(double r[], int ind[], int s, int m);
void hpsort(int n, double *ra, int *ind);


// in mymath4.cpp
//void Jacobi(ComplexMatrix &Ain,double *evals, ComplexMatrix &evecs);

#endif // MYMATH_H


