#ifndef FFT_H
#define FFT_H

#include <complex>
#include <string>

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

namespace ModulePW
{

class FFT
{
public:

	FFT();
	~FFT();
	void initfft(int nx_in, int bigny_in, int nz_in, int lix_in, int rix_in, int ns_in, int nplane_in, 
				 int nproc_in, bool gamma_only_in, bool mpifft_in = false);
	void setupFFT();
	void cleanFFT();

	void fftzfor(std::complex<double>* & in, std::complex<double>* & out);
	void fftzbac(std::complex<double>* & in, std::complex<double>* & out);
	void fftxyfor(std::complex<double>* & in, std::complex<double>* & out);
	void fftxybac(std::complex<double>* & in, std::complex<double>* & out);
	void fftxyr2c(double * &in, std::complex<double>* & out);
	void fftxyc2r(std::complex<double>* & in, double* & out);

#ifdef __MIX_PRECISION
	void executefftwf(std::string instr);
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
	int nxy;
	int bigny;
	int bignxy;
	int lix,rix;// lix: the left edge of the pw ball; rix: the right edge of the pw ball
	int ns; //number of sticks
	int nplane; //number of x-y planes
	int maxgrids; // max between nz * ns and bignxy * nplane
	int nproc; // number of proc.
	std::complex<double> *aux1, *aux2; //fft space, [maxgrids]
	double *r_rspace; //real number space for r, [nplane * nx *ny]
#ifdef __MIX_PRECISION
	std::complex<float> * cf_gspace; //complex number space for g, [ns * nz]
	std::complex<float> * cf_rspace; //complex number space for r, [nplane * nx *ny]
	float *rf_rspace; //real number space for r, [nplane * nx *ny]
#endif


private:
	bool gamma_only;
	bool destroyp;
	bool mpifft; // if use mpi fft, only used when define __FFTW3_MPI
	// fftw_plan plan2r2c;
	// fftw_plan plan2c2r;
	// fftw_plan plan1for;
	// fftw_plan plan1bac;
	// fftw_plan plan2for;
	// fftw_plan plan2bac;
	fftw_plan planzfor;
	fftw_plan planzbac;
	fftw_plan planxfor;
	fftw_plan planxbac;
	fftw_plan planyfor;
	fftw_plan planybac;
	fftw_plan planyr2c;
	fftw_plan planyc2r;
#ifdef __MIX_PRECISION
	bool destroypf;
	fftwf_plan planf2r2c;
	fftwf_plan planf2c2r;
	fftwf_plan planf1for;
	fftwf_plan planf1bac;
	fftwf_plan planf2for;
	fftwf_plan planf2bac;
#endif

};
}

#endif

