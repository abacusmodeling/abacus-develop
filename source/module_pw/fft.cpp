#include "fft.h"
namespace ModulePW
{

FFT::FFT()
{
}

FFT::~FFT()
{
	this->clear();
}
void FFT::clear()
{
	this->cleanFFT();
	if(aux1!=nullptr) {fftw_free(aux1); aux1 = nullptr;}
	if(aux2!=nullptr) {fftw_free(aux2); aux2 = nullptr;}
	if(r_rspace!=nullptr) {fftw_free(r_rspace); r_rspace = nullptr;}
#ifdef __MIX_PRECISION
	this->cleanfFFT();
	if(auxf1!=nullptr) {fftw_free(auxf1); auxf1 = nullptr;}
	if(auxf2!=nullptr) {fftw_free(auxf2); auxf2 = nullptr;}
	if(rf_rspace!=nullptr) {fftw_free(rf_rspace); rf_rspace = nullptr;}
#endif
}

void FFT:: initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, int nproc_in, bool gamma_only_in, bool xprime_in, bool mpifft_in)
{
	this->gamma_only = gamma_only_in;
	this->xprime = xprime_in;
	this->fftnx = this->nx = nx_in;
	this->fftny = this->ny = ny_in;
	if(this->gamma_only)
	{
		if(xprime) 	this->fftnx = int(nx/2) +1;
		else		this->fftny = int(ny/2) +1;
	} 
	this->nz = nz_in;
	this->ns = ns_in;
	this->lixy = lixy_in;
	this->rixy = rixy_in;
	this->nplane = nplane_in;
	this->nproc = nproc_in;
	this->mpifft = mpifft_in;
	this->nxy = this->nx * this-> ny;
	this->fftnxy = this->fftnx * this->fftny;
	this->maxgrids = (this->nz * this->ns > this->nxy * nplane) ? this->nz * this->ns : this->nxy * nplane;
	if(!this->mpifft)
	{
		aux1  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		aux2  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		r_rspace = (double *) fftw_malloc(sizeof(double) * this->nxy * nplane);
#ifdef __MIX_PRECISION
		auxf1  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
		auxf2  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
		rf_rspace = (float *) fftw_malloc(sizeof(float) * this->nxy * nplane);
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
	//                              1 D - Z
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

	//---------------------------------------------------------
	//                              2 D - XY
	//---------------------------------------------------------
	//1D+1D is much faster than 2D FFT!
	
	int *embed = nullptr;
	int npy = this->nplane * this->ny;
	if(this->xprime)
	{
		//It is better to use in-place for stride > 1
		this->planyfor  = fftw_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftw_complex *)aux2, 		  embed, nplane,     1,
				(fftw_complex *)aux2, 	 embed, nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planybac  = fftw_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftw_complex *)aux2, 		  embed, nplane,     1,
				(fftw_complex *)aux2, 	 embed, nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planxr2c  = fftw_plan_many_dft_r2c(  1, &this->nx, npy,	r_rspace , embed, npy,      1,
				(fftw_complex*)aux1, embed, npy,		1,	FFTW_MEASURE   );
			this->planxc2r  = fftw_plan_many_dft_c2r(  1, &this->nx, npy,	(fftw_complex*)aux1 , embed, npy,      1,
				r_rspace, embed, npy,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, npy,     1,
				(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, npy,     1,
				(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
		
	}
	else
	{
		this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftw_complex *)aux2, 		  embed, npy,     1,
				(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftw_complex *)aux2, 		  embed, npy,     1,
				(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->ny, this->nplane,	r_rspace , embed, this->nplane,      1,
				(fftw_complex*)aux1, embed, this->nplane,		1,	FFTW_MEASURE   );
			this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->ny, this->nplane,	(fftw_complex*)aux1 , embed, this->nplane,      1,
				r_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		}
		else
		{

			this->planxfor2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (ny - rixy),	 (fftw_complex *)aux2, 		  embed, npy,     1,
					(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planxbac2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (ny - rixy),	 (fftw_complex *)aux2, 		  embed, npy,     1,
					(fftw_complex *)aux2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
			this->planyfor  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
				(fftw_complex*)aux2, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planybac  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
				(fftw_complex*)aux2, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
	}

	

	destroyp = false;
}

#ifdef __MIX_PRECISION
void FFT :: initplanf()
{
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------

	//               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	//It is better to use out-of-place fft for stride = 1.
	this->planfzfor = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftwf_complex*) auxf1,  &this->nz,  1,  this->nz,
					    (fftwf_complex*) auxf2,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planfzbac = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftwf_complex*) auxf1,  &this->nz,  1,  this->nz,
						(fftwf_complex*) auxf2,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int *embed = nullptr;
	int npy = this->nplane * this->ny;
	if(this->xprime)
	{
		//It is better to use in-place for stride > 1
		this->planfyfor  = fftwf_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftwf_complex *)auxf2, 		  embed, nplane,     1,
				(fftwf_complex *)auxf2, 	 embed, nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfybac  = fftwf_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftwf_complex *)auxf2, 		  embed, nplane,     1,
				(fftwf_complex *)auxf2, 	 embed, nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planfxr2c  = fftwf_plan_many_dft_r2c(  1, &this->nx, npy,	rf_rspace , embed, npy,      1,
				(fftwf_complex*)auxf1, embed, npy,		1,	FFTW_MEASURE   );
			this->planfxc2r  = fftwf_plan_many_dft_c2r(  1, &this->nx, npy,	(fftwf_complex*)auxf1 , embed, npy,      1,
				rf_rspace, embed, npy,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	npy,	 (fftwf_complex *)auxf2, 		  embed, npy,     1,
				(fftwf_complex *)auxf2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	npy,	 (fftwf_complex *)auxf2, 		  embed, npy,     1,
				(fftwf_complex *)auxf2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
		
	}
	else
	{
		this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftwf_complex *)auxf2, 		  embed, npy,     1,
				(fftwf_complex *)auxf2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftwf_complex *)auxf2, 		  embed, npy,     1,
				(fftwf_complex *)auxf2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planfyr2c  = fftwf_plan_many_dft_r2c(  1, &this->ny, this->nplane,	rf_rspace , embed, this->nplane,      1,
				(fftwf_complex*)auxf1, embed, this->nplane,		1,	FFTW_MEASURE   );
			this->planfyc2r  = fftwf_plan_many_dft_c2r(  1, &this->ny, this->nplane,	(fftwf_complex*)auxf1 , embed, this->nplane,      1,
				rf_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planfxfor2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (this->ny - rixy),	 (fftwf_complex *)auxf2, 		  embed, npy,     1,
				(fftwf_complex *)auxf2, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfxbac2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (this->ny - rixy),	 (fftwf_complex *)auxf2, 		  embed, npy,     1,
				(fftwf_complex *)auxf2, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
			this->planfyfor  = fftwf_plan_many_dft(  1, &this->ny, this->nplane,	(fftwf_complex*)auxf2 , embed, this->nplane,      1,
				(fftwf_complex*)auxf2, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfybac  = fftwf_plan_many_dft(  1, &this->ny, this->nplane,	(fftwf_complex*)auxf2 , embed, this->nplane,      1,
				(fftwf_complex*)auxf2, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
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
	if(this->xprime)
	{
		fftw_destroy_plan(planyfor);
		fftw_destroy_plan(planybac);
		if(this->gamma_only)
		{
			fftw_destroy_plan(planxr2c);
			fftw_destroy_plan(planxc2r);
		}
		else
		{
			fftw_destroy_plan(planxfor1);
			fftw_destroy_plan(planxbac1);
		}
	}
	else
	{
		fftw_destroy_plan(planxfor1);
		fftw_destroy_plan(planxbac1);
		if(this->gamma_only)
		{
			fftw_destroy_plan(planyr2c);
			fftw_destroy_plan(planyc2r);
		}
		else
		{
			fftw_destroy_plan(planxfor2);
			fftw_destroy_plan(planxbac2);
			fftw_destroy_plan(planyfor);
			fftw_destroy_plan(planybac);
		}
	}
	destroyp = true;
	return;
}
#ifdef __MIX_PRECISION
void FFT:: cleanfFFT()
{
	if(destroypf==true) return;
	fftwf_destroy_plan(planfzfor);
	fftwf_destroy_plan(planfzbac);
	if(this->xprime)
	{
		fftwf_destroy_plan(planfyfor);
		fftwf_destroy_plan(planfybac);
		if(this->gamma_only)
		{
			fftwf_destroy_plan(planfxr2c);
			fftwf_destroy_plan(planfxc2r);
		}
		else
		{
			fftwf_destroy_plan(planfxfor1);
			fftwf_destroy_plan(planfxbac1);
		}
	}
	else
	{
		fftwf_destroy_plan(planfxfor1);
		fftwf_destroy_plan(planfxbac1);
		if(this->gamma_only)
		{
			fftwf_destroy_plan(planfyr2c);
			fftwf_destroy_plan(planfyc2r);
		}
		else
		{
			fftwf_destroy_plan(planfxfor2);
			fftwf_destroy_plan(planfxbac2);
			fftwf_destroy_plan(planfyfor);
			fftwf_destroy_plan(planfybac);
		}
	}
	destroypf = true;
	return;
}
#endif


void FFT::fftzfor(std::complex<double>* & in, std::complex<double>* & out)
{
	fftw_execute_dft(this->planzfor,(fftw_complex *)in,(fftw_complex *)out);
	return;
}

void FFT::fftzbac(std::complex<double>* & in, std::complex<double>* & out)
{
	fftw_execute_dft(this->planzbac,(fftw_complex *)in, (fftw_complex *)out);
	return;
}

void FFT::fftxyfor(std::complex<double>* & in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftw_execute_dft( this->planxfor1, (fftw_complex *)in, (fftw_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*npy], (fftw_complex *)&out[i*npy]);
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*npy], (fftw_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*npy], (fftw_complex *)&out[i*npy]);
		}

		fftw_execute_dft( this->planxfor1, (fftw_complex *)in, (fftw_complex *)out);
		fftw_execute_dft( this->planxfor2, (fftw_complex *)&in[rixy*nplane], (fftw_complex *)&out[rixy*nplane]);
	}
	return;
}

void FFT::fftxybac(std::complex<double>* & in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
		}

		fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)out);
	}
	else
	{
		fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)out);
		fftw_execute_dft( this->planxbac2, (fftw_complex *)&in[rixy*nplane], (fftw_complex *)&out[rixy*nplane]);

		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
		}
	}
	return;
}

void FFT::fftxyr2c(double* &in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftw_execute_dft_r2c( this->planxr2c, in, (fftw_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&out[i*npy], (fftw_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft_r2c( this->planyr2c, &in[i*npy], (fftw_complex*)&out[i*npy] );
		}

		fftw_execute_dft( this->planxfor1, (fftw_complex *)out, (fftw_complex *)out);
	}
	return;
}


void FFT::fftxyc2r(std::complex<double>* & in, double* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex *)&in[i*npy], (fftw_complex *)&in[i*npy]);
		}

		fftw_execute_dft_c2r( this->planxc2r, (fftw_complex *)in, out);
	}
	else
	{
		fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)in);

		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft_c2r( this->planyc2r, (fftw_complex*)&in[i*npy], &out[i*npy] );
		}
	}
	return;
}


#ifdef __MIX_PRECISION
void FFT::fftfzfor(std::complex<float>* & in, std::complex<float>* & out)
{
	fftwf_execute_dft(this->planfzfor,(fftwf_complex *)in,(fftwf_complex *)out);
	return;
}

void FFT::fftfzbac(std::complex<float>* & in, std::complex<float>* & out)
{
	fftwf_execute_dft(this->planfzbac,(fftwf_complex *)in, (fftwf_complex *)out);
	return;
}

void FFT::fftfxyfor(std::complex<float>* & in, std::complex<float>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)in, (fftwf_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&out[i*npy]);
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&out[i*npy]);
		}

		fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)in, (fftwf_complex *)out);
		fftwf_execute_dft( this->planfxfor2, (fftwf_complex *)&in[rixy*nplane], (fftwf_complex *)&out[rixy*nplane]);
	}
	return;
}

void FFT::fftfxybac(std::complex<float>* & in, std::complex<float>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*npy], (fftwf_complex*)&out[i*npy] );
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*npy], (fftwf_complex*)&out[i*npy] );
		}

		fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)out);
	}
	else
	{
		fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)out);
		fftwf_execute_dft( this->planfxbac2, (fftwf_complex *)&in[rixy*nplane], (fftwf_complex *)&out[rixy*nplane]);

		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*npy], (fftwf_complex*)&out[i*npy] );
		}
	}
	return;
}

void FFT::fftfxyr2c(float* &in, std::complex<float>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftwf_execute_dft_r2c( this->planfxr2c, in, (fftwf_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&out[i*npy], (fftwf_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft_r2c( this->planfyr2c, &in[i*npy], (fftwf_complex*)&out[i*npy] );
		}

		fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)out, (fftwf_complex *)out);
	}
	
	return;
}


void FFT::fftfxyc2r(std::complex<float>* & in, float* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&in[i*npy]);
		}

		fftwf_execute_dft_c2r( this->planfxc2r, (fftwf_complex *)in, out);
	}
	else
	{
		fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)in);

		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft_c2r( this->planfyc2r, (fftwf_complex*)&in[i*npy], &out[i*npy] );
		}
	}
	return;
}
#endif
}
