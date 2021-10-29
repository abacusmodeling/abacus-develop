#include "fft.h"

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
	
void FFT:: setupFFT3D()
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

	void FFT:: cleanFFT3D()
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

		destroyp == true;
		return;
	}

