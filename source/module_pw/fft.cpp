#include "fft.h"
#include "../module_base/tool_quit.h"
namespace ModulePW
{

FFT::FFT()
{
	bignx = nx = ny = nz = 2;
	bignxy = nxy = 4;
	ns = 1; 
	nplane = 1;
	mpifft = false; 
	destroyp = true;
	c_rspace = c_gspace = NULL;
	r_rspace = NULL;
#ifdef __MIX_PRECISION
	destroypf = true;
	cf_rspace = cf_gspace = NULL;
	rf_rspace = NULL;
#endif
}

FFT::~FFT()
{
	this->cleanFFT();
	if(c_gspace!=NULL) fftw_free(c_gspace);
	if(c_rspace!=NULL) fftw_free(c_rspace);
	//if(r_gspace!=NULL) fftw_free(r_gspace);
	if(r_rspace!=NULL) fftw_free(r_rspace);
#ifdef __MIX_PRECISION
	if(cf_gspace!=NULL) fftw_free(cf_gspace);
	if(cf_rspace!=NULL) fftw_free(cf_rspace);
	//if(rf_gspace!=NULL) fftw_free(rf_gspace);
	if(rf_rspace!=NULL) fftw_free(rf_rspace);
#endif
}

void FFT:: initfft(int bignx_in, int nx_in, int ny_in , int nz_in, int ns_in, int nplane_in, bool mpifft_in)
{
	this->bignx = bignx_in;
	this->nx = nx_in;
	this->ny = ny_in;
	this->nz = nz_in;
	this->ns = ns_in;
	this->nplane = nplane_in;
	this->mpifft = mpifft_in;
	this->nxy = this->nx * this-> ny;
	this->bignxy = this->bignx * this-> ny;
	if(!this->mpifft)
	{
		c_gspace  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * this->nz * this->ns);
		//r_gspace = (double *) fftw_malloc(sizeof(double) * this->nz * this->ns);
		c_rspace  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * this->bignxy * nplane);
		r_rspace = (double *) fftw_malloc(sizeof(double) * this->bignxy * nplane);
#ifdef __MIX_PRECISION
		cf_gspace  = (std::complex<float> *)fftw_malloc(sizeof(fftwf_complex) * this->nz * this->ns);
		//rf_gspace = (float *)fftw_malloc(sizeof(float) * this->nz * this->ns);
		cf_rspace  = (std::complex<float> *)fftw_malloc(sizeof(fftwf_complex) * this->bignxy * nplane);
		rf_rspace = (float *)fftw_malloc(sizeof(float) * this->bignxy * nplane);
#endif
	}
	else
	{
		
	}
	
}

void FFT:: setupFFT()
{
	if(!this->mpifft)
	{
		this->initplan();
#ifdef __MIX_PRECISION
		this->initplanf();
#endif
	}
#if defined(__FFTW3_MPI) && defined(__MPI)
	else
	{
		this->initplan_mpi();
#ifdef __MIX_PRECISION
		this->initplanf_mpi();
#endif
	}
#endif
	return;
}
	
void FFT :: initplan()
{
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------

	//               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	this->plan1for = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftw_complex*) c_gspace,  &this->nz,  1,  this->nz,
					    (fftw_complex*) c_gspace,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->plan1bac = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftw_complex*) c_gspace,  &this->nz,  1,  this->nz,
						(fftw_complex*) c_gspace,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	// this->plan1r2c = fftw_plan_many_dft_r2c( 1,    &this->nz,  this->ns,  
	// 									r_gspace,  &this->nz,  1,  this->nz,
	// 					(fftw_complex*) c_gspace,  &this->nz,  1,  this->nz,  FFTW_MEASURE);
	
	// this->plan1c2r = fftw_plan_many_dft_c2r( 1,   &this->nz,  this->ns,  
	// 					(fftw_complex*) c_gspace, &this->nz,  1,  this->nz,
	// 									r_gspace, &this->nz,  1,  this->nz,  FFTW_MEASURE);

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int nrank[2] = {this->ny,this->bignx};
	this->plan2for = fftw_plan_many_dft(       2,   nrank,  this->nplane,  
						(fftw_complex*) c_rspace,   nrank,  this->nplane,   1,
						(fftw_complex*) c_rspace,   nrank,  this->nplane,   1,  FFTW_FORWARD,  FFTW_MEASURE);

	this->plan2bac = fftw_plan_many_dft(       2,   nrank,  this->nplane,  
						(fftw_complex*) c_rspace,   nrank,  this->nplane,   1,
						(fftw_complex*) c_rspace,   nrank,  this->nplane,   1,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	this->plan2r2c = fftw_plan_many_dft_r2c(   2,   nrank,  this->nplane,  
										r_rspace,   nrank,  this->nplane,   1,
						(fftw_complex*) c_rspace,   nrank,  this->nplane,   1,  FFTW_MEASURE);
	
	this->plan2c2r = fftw_plan_many_dft_c2r(   2,   nrank,  this->nplane,  
						(fftw_complex*) c_rspace,   nrank,  this->nplane,   1,
										r_rspace,   nrank,  this->nplane,   1,  FFTW_MEASURE);
	destroyp = false;
}

#ifdef __MIX_PRECISION
void FFT :: initplanf()
{
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------

	//              fftwf_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	this->planf1for = fftwf_plan_many_dft(     1,  &this->nz,  this->ns,  
						(fftwf_complex*)c_gspace,  &this->nz,  1,  this->nz,
						(fftwf_complex*)c_gspace,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planf1bac = fftwf_plan_many_dft(     1,  &this->nz,  this->ns,  
						(fftwf_complex*)c_gspace,  &this->nz,  1,  this->nz,
						(fftwf_complex*)c_gspace,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	// this->planf1r2c = fftwf_plan_many_dft_r2c( 1,  &this->nz,  this->ns,  
	// 									r_gspace,  &this->nz,  1,  this->nz,
	// 					(fftwf_complex*)c_gspace,  &this->nz,  1,  this->nz,  FFTW_MEASURE);
	
	// this->planf1c2r = fftwf_plan_many_dft_c2r( 1, &this->nz,  this->ns,  
	// 					(fftwf_complex*)c_gspace, &this->nz,  1,  this->nz,
	// 									r_gspace, &this->nz,  1,  this->nz,  FFTW_MEASURE);

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int nrank[2] = {this->ny,this->bignx};
	this->planf2for = fftwf_plan_many_dft(     2,   nrank,  this->nplane,  
						(fftwf_complex*)c_rspace,   nrank,  this->nplane,   1,
						(fftwf_complex*)c_rspace,   nrank,  this->nplane,   1,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planf2bac = fftwf_plan_many_dft(     2,   nrank,  this->nplane,  
						(fftwf_complex*)c_rspace,   nrank,  this->nplane,   1,
						(fftwf_complex*)c_rspace,   nrank,  this->nplane,   1,  FFTW_BACKWARD,  FFTW_MEASURE);
	
	this->planf2r2c = fftwf_plan_many_dft_r2c( 2,   nrank,  this->nplane,  
										r_rspace,   nrank,  this->nplane,   1,
						(fftwf_complex*)c_rspace,   nrank,  this->nplane,   1,  FFTW_MEASURE);
	
	this->planf2c2r = fftwf_plan_many_dft_c2r( 2,   nrank,  this->nplane,  
						(fftwf_complex*)c_rspace,   nrank,  this->nplane,   1,
										r_rspace,   nrank,  this->nplane,   1,  FFTW_MEASURE);
	destroypf = false;
}
#endif

void FFT :: initplan_mpi()
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
	// fftw_destroy_plan(plan1r2c);
	// fftw_destroy_plan(plan1c2r);
	fftw_destroy_plan(plan2for);
	fftw_destroy_plan(plan2bac);
	fftw_destroy_plan(plan2r2c);
	fftw_destroy_plan(plan2c2r);
	destroyp == true;

#ifdef __MIX_PRECISION
	if(destroypf==true) return;
	fftw_destroy_plan(planf1for);
	fftw_destroy_plan(planf1bac);
	// fftw_destroy_plan(planf1r2c);
	// fftw_destroy_plan(planf1c2r);
	fftw_destroy_plan(planf2for);
	fftw_destroy_plan(planf2bac);
	fftw_destroy_plan(planf2r2c);
	fftw_destroy_plan(planf2c2r);
	destroypf == true;
#endif

	return;
}

void FFT::executefftw(std::string instr)
{
	if(instr == "1for")
		fftw_execute(this->plan1for);
	else if(instr == "2for")
		fftw_execute(this->plan2for);
	else if(instr == "1bac")
		fftw_execute(this->plan1bac);
	else if(instr == "2bac")
		fftw_execute(this->plan2bac);
	// else if(instr == "1r2c")
	// 	fftw_execute(this->plan1r2c);
	else if(instr == "2r2c")
		fftw_execute(this->plan2r2c);
	// else if(instr == "1c2r")
	// 	fftw_execute_dft(this->plan1c2r);
	else if(instr == "2c2r")
		fftw_execute(this->plan2c2r);
	else
	{
		ModuleBase::WARNING_QUIT("FFT", "Wrong input for excutefftw");
	}
}

#ifdef __MIX_PRECISION
void executefftwf(std::string instr)
{
	if(instr == "1for")
		fftwf_execute(this->planf1for);
	else if(instr == "2for")
		fftwf_execute(this->planf2for);
	else if(instr == "1bac")
		fftwf_execute(this->planf1bac);
	else if(instr == "2bac")
		fftwf_execute(this->planf2bac);
	// else if(instr == "1r2c")
	// 	fftwf_execute(this->planf1r2c);
	else if(instr == "2r2c")
		fftwf_execute(this->planf2r2c);
	// else if(instr == "1c2r")
	// 	fftwf_execute(this->planf1c2r);
	else if(instr == "2c2r")
		fftwf_execute(this->planf2c2r);
	else
	{
		ModuleBase::WARNING_QUIT("FFT", "Wrong input for excutefftwf");
	}
}
#endif
}
