#include "fft.h"
#include <iostream> 

namespace ModulePW
{

FFT::FFT()
{
	bigny = nx = ny = nz = 2;
	bignxy = nxy = 4;
	ns = 1; 
	nplane = 1;
	mpifft = false; 
	destroyp = true;
	gamma_only = false;
	aux2 = aux1 = NULL;
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
	if(aux1!=NULL) fftw_free(aux1);
	if(aux2!=NULL) fftw_free(aux2);
#ifdef __MIX_PRECISION
	if(cf_gspace!=NULL) fftw_free(cf_gspace);
	if(cf_rspace!=NULL) fftw_free(cf_rspace);
	if(rf_rspace!=NULL) fftw_free(rf_rspace);
#endif
}

void FFT:: initfft(int nx_in, int bigny_in, int nz_in, int lix_in, int rix_in, int ns_in, int nplane_in, int nproc_in, bool gamma_only_in, bool mpifft_in)
{
	this->gamma_only = gamma_only_in;
	this->nx = nx_in;
	this->bigny = bigny_in;
	if(this->gamma_only) this->ny = int(bigny/2) +1;
	else	this->ny = this->bigny;
	this->nz = nz_in;
	this->ns = ns_in;
	this->lix = lix_in;
	this->rix = rix_in;
	this->nplane = nplane_in;
	this->nproc = nproc_in;
	this->mpifft = mpifft_in;
	this->nxy = this->nx * this-> ny;
	this->bignxy = this->nx * this->bigny;
	this->maxgrids = (this->nz * this->ns > this->bignxy * nplane) ? this->nz * this->ns : this->bignxy * nplane;
	if(!this->mpifft)
	{
		aux1  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		aux2  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		r_rspace = (double *) aux1;
#ifdef __MIX_PRECISION
		cf_gspace  = (std::complex<float> *)fftw_malloc(sizeof(fftwf_complex) * this->nz * this->ns);
		cf_rspace  = (std::complex<float> *)fftw_malloc(sizeof(fftwf_complex) * this->bignxy * nplane);
		if(this->gamma_only)
		{
			rf_rspace = (float *)fftw_malloc(sizeof(float) * this->bignxy * nplane);
		}
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
	
	//It is better to use out-of-place fft for stride = 1.
	this->planzfor = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftw_complex*) aux1,  &this->nz,  1,  this->nz,
					    (fftw_complex*) aux2,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planzbac = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftw_complex*) aux1,  &this->nz,  1,  this->nz,
						(fftw_complex*) aux2,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);

	// this->planzfor = fftw_plan_dft_1d(this->nz,(fftw_complex*) aux1,(fftw_complex*) aux1, FFTW_FORWARD,  FFTW_MEASURE);
	// this->planzbac = fftw_plan_dft_1d(this->nz,(fftw_complex*) aux1,(fftw_complex*) aux1,FFTW_BACKWARD,  FFTW_MEASURE);

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	//int nrank[2] = {this->nx,this->bigny};
	int *embed = NULL;
	// It seems 1D+1D is much faster than 2D FFT!
	if(this->gamma_only)
	{
		// int padnpy = this->nplane * this->ny * 2;
		// int rankc[2] = {this->nx, this->padnpy};
		// int rankd[2] = {this->nx, this->padnpy*2};
		// // It seems 1D+1D is much faster than 2D FFT!
		// this->plan2r2c = fftw_plan_many_dft_r2c(   2,   nrank,  this->nplane,  
		// 									r_rspace,   rankd,  this->nplane,   1,
		// 					(fftw_complex*) aux2,   rankc,  this->nplane,   1,  FFTW_MEASURE);

		// this->plan2c2r = fftw_plan_many_dft_c2r(   2,   nrank,  this->nplane,  
		// 					(fftw_complex*) aux2,   rankc,  this->nplane,   1,
		// 									r_rspace,   rankd,  this->nplane,   1,  FFTW_MEASURE);

		int npy = this->nplane * this->ny;
		int bignpy = this->nplane * this->bigny;
		this->planxfor  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->bigny, this->nplane,	r_rspace , embed, this->nplane,      1,
			(fftw_complex*)aux1, embed, this->nplane,		1,	FFTW_MEASURE   );
		this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->bigny, this->nplane,	(fftw_complex*)aux1 , embed, this->nplane,      1,
			r_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		
		// int padnpy = npy * 2;
		// this->planxfor  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, padnpy,     1,
		// 		(fftw_complex *)aux2, 	 embed, padnpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		// this->planxbac  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, padnpy,     1,
		// 		(fftw_complex *)aux2, 	 embed, padnpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		// this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->bigny, this->nplane,	r_rspace , embed, this->nplane*2,      1,
		// 	(fftw_complex*)aux2, embed, this->nplane,		1,	FFTW_MEASURE   );
		// this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->bigny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
		// 	r_rspace, embed, this->nplane*2,		1,			FFTW_MEASURE   );

	}
	else
	{
		// 	this->plan2for = fftw_plan_many_dft(       2,   nrank,  this->nplane,  
		// 						(fftw_complex*) aux2,   embed,  this->nplane,   1,
		// 						(fftw_complex*) aux2,  embed,  this->nplane,   1,  FFTW_FORWARD,  FFTW_MEASURE);

		// 	this->plan2bac = fftw_plan_many_dft(       2,   nrank,  this->nplane,  
		// 						(fftw_complex*) aux2,   embed,  this->nplane,   1,
		// 						(fftw_complex*) aux2,   embed,  this->nplane,   1,  FFTW_BACKWARD,  FFTW_MEASURE);

		//It is better to use in-place for stride > 1
		int npy = this->nplane * this->ny;
		this->planxfor  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, npy,     1,
				(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, npy,     1,
				(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		this->planyfor  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
			(fftw_complex*)aux2, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planybac  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
			(fftw_complex*)aux2, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
	}

	

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
						(fftwf_complex*)aux1,  &this->nz,  1,  this->nz,
						(fftwf_complex*)aux1,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planf1bac = fftwf_plan_many_dft(     1,  &this->nz,  this->ns,  
						(fftwf_complex*)aux1,  &this->nz,  1,  this->nz,
						(fftwf_complex*)aux1,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int nrank[2] = {this->nx,this->bigny};
	
	if(this->gamma_only)
	{
		this->planf2r2c = fftwf_plan_many_dft_r2c( 2,   nrank,  this->nplane,  
											r_rspace,   nrank,  this->nplane,   1,
							(fftwf_complex*)aux2,   nrank,  this->nplane,   1,  FFTW_MEASURE);

		this->planf2c2r = fftwf_plan_many_dft_c2r( 2,   nrank,  this->nplane,  
							(fftwf_complex*)aux2,   nrank,  this->nplane,   1,
											r_rspace,   nrank,  this->nplane,   1,  FFTW_MEASURE);
	}
	else
	{
		this->planf2for = fftwf_plan_many_dft(     2,   nrank,  this->nplane,  
							(fftwf_complex*)aux2,   nrank,  this->nplane,   1,
							(fftwf_complex*)aux2,   nrank,  this->nplane,   1,  FFTW_FORWARD,  FFTW_MEASURE);

		this->planf2bac = fftwf_plan_many_dft(     2,   nrank,  this->nplane,  
							(fftwf_complex*)aux2,   nrank,  this->nplane,   1,
							(fftwf_complex*)aux2,   nrank,  this->nplane,   1,  FFTW_BACKWARD,  FFTW_MEASURE);
	}
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
	fftw_destroy_plan(planzfor);
	fftw_destroy_plan(planzbac);
	fftw_destroy_plan(planxfor);
	fftw_destroy_plan(planxbac);
	if(this->gamma_only)
	{
		fftw_destroy_plan(planyr2c);
		fftw_destroy_plan(planyc2r);
	}
	else
	{
		fftw_destroy_plan(planyfor);
		fftw_destroy_plan(planybac);
	}
	destroyp = true;

#ifdef __MIX_PRECISION
	if(destroypf==true) return;
	fftw_destroy_plan(planf1for);
	fftw_destroy_plan(planf1bac);
	if(this->gamma_only)
	{
		fftw_destroy_plan(planf2r2c);
		fftw_destroy_plan(planf2c2r);
	}
	else
	{
		fftw_destroy_plan(planf2for);
		fftw_destroy_plan(planf2bac);
	}
	destroypf = true;
#endif

	return;
}

void FFT::fftzfor(std::complex<double>* & in, std::complex<double>* & out)
{
	// for(int i = 0 ; i < this->ns ; ++i)
	// {
	// 	fftw_execute_dft(this->planzfor,(fftw_complex *)&in[i*nz],(fftw_complex *)&out[i*nz]);
	// }
	fftw_execute_dft(this->planzfor,(fftw_complex *)in,(fftw_complex *)out);
	return;
}

void FFT::fftzbac(std::complex<double>* & in, std::complex<double>* & out)
{
	// for(int i = 0 ; i < this->ns ; ++i)
	// {
	// 	fftw_execute_dft(this->planzbac,(fftw_complex *)&aux1[i*nz],(fftw_complex *)&aux1[i*nz]);
	// }
	fftw_execute_dft(this->planzbac,(fftw_complex *)in, (fftw_complex *)out);
	return;
}

void FFT::fftxyfor(std::complex<double>* & in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	fftw_execute_dft( this->planxfor, (fftw_complex *)in, (fftw_complex *)out);
	for (int i=0; i<=this->lix;++i)
	{
		fftw_execute_dft( this->planyfor, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
	}
	for (int i=this->rix; i<this->nx;++i)
	{
		fftw_execute_dft( this->planyfor, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
	}
	return;
}

void FFT::fftxybac(std::complex<double>* & in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
		
	for (int i=0; i<=this->lix;++i)
	{
		fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
	}
	for (int i=this->rix; i<this->nx;++i)
	{
		fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
	}
	fftw_execute_dft( this->planxbac, (fftw_complex *)in, (fftw_complex *)out);
	return;
}

void FFT::fftxyr2c(double* &in, std::complex<double>* & out)
{
	//int npy = this->nplane * this-> ny;
	int bignpy = this->nplane * this-> bigny;
	// int padnpy = this->nplane * this-> ny * 2;
	for (int i=0; i<this->nx;++i)
	{
		fftw_execute_dft_r2c( this->planyr2c, &in[i*bignpy*2], (fftw_complex*)&out[i*bignpy] );
		if((double*)&out[i*bignpy] != &in[i*bignpy*2] ) std::cout<<i<<"wrond\n";
		// fftw_execute_dft_r2c( this->planyfor, &r_rspace[4*i*padnpy], (fftw_complex*)&aux2[i*padnpy] );
	}
	fftw_execute_dft( this->planxfor, (fftw_complex *)out, (fftw_complex *)out);
	return;
}


void FFT::fftxyc2r(std::complex<double>* & in, double* & out)
{
	//int npy = this->nplane * this-> ny;
	int bignpy = this->nplane * this-> bigny;
	// int padnpy = this->nplane * this-> ny * 2;
	fftw_execute_dft( this->planxbac, (fftw_complex *)in, (fftw_complex *)in);
	for (int i=0; i<this->nx;++i)
	{
		fftw_execute_dft_c2r( this->planyc2r, (fftw_complex*)&in[i*bignpy], &out[i*bignpy*2] );
		// fftw_execute_dft_c2r( this->planybac, (fftw_complex*)&aux2[i*padnpy], &r_rspace[4*i*padnpy] );
	}
	return;
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
	else if(instr == "2r2c")
		fftwf_execute(this->planf2r2c);
	else if(instr == "2c2r")
		fftwf_execute(this->planf2c2r);
}
#endif
}
