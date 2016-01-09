/*******************************************************
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
*********************************************************/

/* mymath.h file */

#ifndef MYMATH_H
#define MYMATH_H
#include "../src_pw/tools.h"

#ifdef __FFTW2
#include "../src_parallel/fftw.h"
#elif defined __FFTW3
#include "../src_parallel/fftw3.h"
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


