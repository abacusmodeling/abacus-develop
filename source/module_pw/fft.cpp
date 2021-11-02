#include "fft.h"
#include "../module_base/tool_quit.h"

FFT::FFT()
{
	nx = ny = nz = 2;
	int ns = 1; 
	int nplane = 1;
	int ffttype = 1; 
	destroyp = true;
}

FFT::~FFT()
{
	this->cleanFFT3D();
}

void FFT:: initfft(int nx_in, int ny_in , int nz_in, int ns_in, int nplane_in, int ffttype_in)
{
	this->nx = nx_in;
	this->ny = ny_in;
	this->nz = nz_in;
	this->ns = ns_in;
	this->nplane = nplane_in;
	this->ffttype = ffttype_in;
}

void FFT:: setupSFFT()
{
	if(this->ffttype == 1)
	{
		this->initpland();
#ifdef __MIX_PRECISION
		this->initplanf();
#endif
	}
#if defined(__FFTW3_MPI) && defined(__MPI)
	else if(this->ffttype == 2)
	{
		this->initpland_mpi();
#ifdef __MIX_PRECISION
		this->initplanf_mpi();
#endif
	}
#endif
	else
	{
		ModuleBase::WARNING_QUIT("FFT", "No such fft type!");
	}
	return;
}
	
void FFT :: initpland()
{
	fftw_complex * tmp_in;
	fftw_complex * tmp_out;
	double *tmp_r;
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------
	tmp_in  = fftw_malloc(sizeof(fftw_complex) * nz * ns);
	tmp_out = fftw_malloc(sizeof(fftw_complex) * nz * ns);
	tmp_r = fftw_malloc(sizeof(double) * nz * ns);

	//               fftw_plan_many_dft(int rank, const int *n, int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	this->plan1for = fftw_plan_many_dft( 1,   &this->nz,  this->ns,  
										tmp_in,   &this->nz,  1,  this->nz,
										tmp_out,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->plan1bac = fftw_plan_many_dft( 1,   &this->nz,  this->ns,  
										tmp_in,   &this->nz,  1,  this->nz,
										tmp_out,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	this->plan1r2c = fftw_plan_many_dft_r2c( 1,   &this->nz,  this->ns,  
										tmp_r,   &this->nz,  1,  this->nz,
										tmp_out,  &this->nz,  1,  this->nz,  FFTW_MEASURE);
	
	this->plan1c2r = fftw_plan_many_dft_c2r( 1,   &this->nz,  this->ns,  
										tmp_in,  &this->nz,  1,  this->nz,
										tmp_r,   &this->nz,  1,  this->nz,  FFTW_MEASURE);
	fftw_free(tmp_in);
	fftw_free(tmp_out);
	fftw_free(tmp_r);

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	int nxy = this->nx * this-> ny;
	tmp_in  = fftw_malloc(sizeof(fftw_complex) * nxy * nplane);
	tmp_out = fftw_malloc(sizeof(fftw_complex) * nxy * nplane);
	tmp_r = fftw_malloc(sizeof(double) * nz * ns);
	
	int * nrank = {2, 2};
	this->plan1for = fftw_plan_many_dft( 2,   nrank,  this->nz,  
										tmp_in,   nrank,  1,  nxy,
										tmp_out,  nrank,  1,  nxy,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->plan1bac = fftw_plan_many_dft( 2,   nrank,  this->nz,  
										tmp_in,   nrank,  1,  nxy,
										tmp_out,  nrank,  1,  nxy,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	this->plan1r2c = fftw_plan_many_dft_r2c( 2,   nrank,  this->nz,  
										tmp_r,    nrank,  1,  nxy,
										tmp_out,  nrank,  1,  nxy,  FFTW_MEASURE);
	
	this->plan1c2r = fftw_plan_many_dft_c2r( 2,   nrank,  this->nz,  
										tmp_in,   nrank,  1,  nxy,
										tmp_r,    nrank,  1,  nxy,  FFTW_MEASURE);

	fftw_free(tmp_in);
	fftw_free(tmp_out);
	fftw_free(tmp_r);


	destroyp = false;
}

#ifdef __MIX_PRECISION
void FFT :: initplanf()
{
	fftwf_complex * tmp_in;
	fftwf_complex * tmp_out;
	float *tmp_r;
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------
	tmp_in  = fftw_malloc(sizeof(fftwf_complex) * nz * ns);
	tmp_out = fftw_malloc(sizeof(fftwf_complex) * nz * ns);
	tmp_r = fftw_malloc(sizeof(float) * nz * ns);

	//               fftwf_plan_many_dft(int rank, const int *n, int howmany,
	//					                fftwf_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftwf_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	this->planf1for = fftwf_plan_many_dft( 1,   &this->nz,  this->ns,  
										tmp_in,   &this->nz,  1,  this->nz,
										tmp_out,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planf1bac = fftwf_plan_many_dft( 1,   &this->nz,  this->ns,  
										tmp_in,   &this->nz,  1,  this->nz,
										tmp_out,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	this->planf1r2c = fftwf_plan_many_dft_r2c( 1,   &this->nz,  this->ns,  
										tmp_r,   &this->nz,  1,  this->nz,
										tmp_out,  &this->nz,  1,  this->nz,  FFTW_MEASURE);
	
	this->planf1c2r = fftwf_plan_many_dft_c2r( 1,   &this->nz,  this->ns,  
										tmp_in,  &this->nz,  1,  this->nz,
										tmp_r,   &this->nz,  1,  this->nz,  FFTW_MEASURE);
	fftw_free(tmp_in);
	fftw_free(tmp_out);
	fftw_free(tmp_r);

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	int nxy = this->nx * this-> ny;
	tmp_in  = fftw_malloc(sizeof(fftwf_complex) * nxy * nplane);
	tmp_out = fftw_malloc(sizeof(fftwf_complex) * nxy * nplane);
	tmp_r = fftw_malloc(sizeof(float) * nz * ns);
	
	int * nrank = {2, 2};
	this->planf1for = fftwf_plan_many_dft( 2,   nrank,  this->nz,  
										tmp_in,   nrank,  1,  nxy,
										tmp_out,  nrank,  1,  nxy,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planf1bac = fftwf_plan_many_dft( 2,   nrank,  this->nz,  
										tmp_in,   nrank,  1,  nxy,
										tmp_out,  nrank,  1,  nxy,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	this->planf1r2c = fftwf_plan_many_dft_r2c( 2,   nrank,  this->nz,  
										tmp_r,    nrank,  1,  nxy,
										tmp_out,  nrank,  1,  nxy,  FFTW_MEASURE);
	
	this->planf1c2r = fftwf_plan_many_dft_c2r( 2,   nrank,  this->nz,  
										tmp_in,   nrank,  1,  nxy,
										tmp_r,    nrank,  1,  nxy,  FFTW_MEASURE);

	fftw_free(tmp_in);
	fftw_free(tmp_out);
	fftw_free(tmp_r);


	destroyp = false;
}
#endif

void FFT :: initpland_mpi()
{

}

#ifdef __MIX_PRECISION
void FFT :: initplanf_mpi()
{
	
}
#endif

void FFT:: cleanFFT()
{
	if(destroyp==true) return;
	fftw_destroy_plan(plan1for);
	fftw_destroy_plan(plan1bac);
	fftw_destroy_plan(plan1r2c);
	fftw_destroy_plan(plan1c2r);
	fftw_destroy_plan(plan2for);
	fftw_destroy_plan(plan2bac);
	fftw_destroy_plan(plan2r2c);
	fftw_destroy_plan(plan2c2r);
#ifdef __MIX_PRECISION
	fftw_destroy_plan(planf1for);
	fftw_destroy_plan(planf1bac);
	fftw_destroy_plan(planf1r2c);
	fftw_destroy_plan(planf1c2r);
	fftw_destroy_plan(planf2for);
	fftw_destroy_plan(planf2bac);
	fftw_destroy_plan(planf2r2c);
	fftw_destroy_plan(planf2c2r);
#endif
	destroyp == true;
	return;
}

void FFT:: executefor(fftw_complex *in, fftw_complex* out, int n)
{
	if(n == 1)
	{
		fftw_execute_dft(this->plan1for, in , out)
	}
	else if(n==2)
	{
		fftw_execute_dft(this->plan2for, in , out)
	}
	else
	{
		ModuleBase::WARNING_QUIT("FFT", "We can only calculate 1D or 2D FFT");
	}
	return;
}


void FFT:: executebac(fftw_complex *in, fftw_complex* out, int n)
{
	if(n == 1)
	{
		fftw_execute_dft(this->plan1bac, in , out)
	}
	else if(n==2)
	{
		fftw_execute_dft(this->plan2bac, in , out)
	}
	else
	{
		ModuleBase::WARNING_QUIT("FFT", "We can only calculate 1D or 2D FFT");
	}
	return;
}

