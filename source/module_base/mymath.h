#ifndef MYMATH_H
#define MYMATH_H
//#include "../src_pw/tools.h"
using namespace std;
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

double rndm();
void simpson(const int mesh,const double *func,const double *rab, double &asum);

double my_erf(double x);
double my_erfc(double x);
double gauss_freq(double x);

void heapsort(int n, double *r, int *ind);
void heapAjust(double r[], int ind[], int s, int m);
void hpsort(int n, double *ra, int *ind);

#endif // MYMATH_H


