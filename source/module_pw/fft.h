#ifndef FFT_H
#define FFT_H

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


class FFT: public Parallel_PW
{
public:

	FFT();
	~FFT();
	void initfft(int nx_in, int ny_in , int nz_in, int ns_in, int nplane_in, int ffttype_in);
	void setupSFFT();
	void initpland();
	void initpland_mpi();
#ifdef __MIX_PRECISION
	void initplanf();
	void initplanf_mpi();
#endif
	void cleanFFT();

	void executefor(fftw_complex *, fftw_complex*, int n);
	void executebac(fftw_complex *, fftw_complex*, int n);
	void executer2c(double *, fftw_complex*, int n);
	void executec2r(fftw_complex *, fftw_complex*, int n);
#ifdef __MIX_PRECISION
	void executeffor(fftwf_complex *, fftwf_complex*, int n);
	void executefbac(fftwf_complex *, fftwf_complex*, int n);
	void executefr2c(float *, fftwf_complex*, int n);
	void executefc2r(fftwf_complex *, float*, int n);
#endif

public:
	int nx,ny,nz;
	int ns; //number of sticks
	int nplane; //number of x-y planes
	int ffttype; // type of FFT


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

