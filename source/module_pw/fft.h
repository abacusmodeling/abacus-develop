#ifndef FFT_H
#define FFT_H

#include <complex>

#include "fftw3.h"
#if defined(__FFTW3_MPI) && defined(__MPI)
#include <fftw3-mpi.h>
//#include "fftw3-mpi_mkl.h"
#endif

#ifdef __MIX_PRECISION
#include "fftw3f.h"
#if defined(__FFTW3_MPI) && defined(__MPI)
#include "fftw3f-mpi.h"
//#include "fftw3-mpi_mkl.h"
#endif
#endif



class FFT
{
public:

	FFT();
	~FFT();
	void initfft(int nx_in, int ny_in , int nz_in, int ns_in, int nplane_in, bool mpifft_in = false);
	void setupFFT();
	void cleanFFT();

	void executefftw(string instr);
#ifdef __MIX_PRECISION
	void executefftwf(string instr);
#endif

private:
	void initplan();
	void initplan_mpi();
#ifdef __MIX_PRECISION
	void initplanf();
	void initplanf_mpi();
#endif
	
public:
	int nx,ny,nz;
	int nxy,nxyz;
	int ns; //number of sticks
	int nplane; //number of x-y planes
	bool mpifft; // if use mpi fft, only used when define __FFTW3_MPI
	complex<double> * c_gspace; //complex number space for g, [ns * nz]
	complex<double> * c_rspace; //complex number space for r, [nplane * nx *ny]
	double *r_gspace; //real number space for g, [ns * nz]
	double *r_rspace; //real number space for r, [nplane * nx *ny]
#ifdef __MIX_PRECISION
	complex<float> * cf_gspace; //complex number space for g, [ns * nz]
	complex<float> * cf_rspace; //complex number space for r, [nplane * nx *ny]
	float *rf_gspace; //real number space for g, [ns * nz]
	float *rf_rspace; //real number space for r, [nplane * nx *ny]
#endif


private:
	bool destroyp;
	fftw_plan plan1_r2c;
	fftw_plan plan1_c2r;
	fftw_plan plan2_r2c;
	fftw_plan plan2_c2r;
	fftw_plan plan1_for;
	fftw_plan plan1_bac;
	fftw_plan plan2_for;
	fftw_plan plan2_bac;
#ifdef __MIX_PRECISION
	bool destroypf;
	fftwf_plan planf1_r2c;
	fftwf_plan planf1_c2r;
	fftwf_plan planf2_r2c;
	fftwf_plan planf2_c2r;
	fftwf_plan planf1_for;
	fftwf_plan planf1_bac;
	fftwf_plan planf2_for;
	fftwf_plan planf2_bac;
#endif

};

#endif

