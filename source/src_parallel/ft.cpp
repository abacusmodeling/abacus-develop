#include "ft.h"
#include "parallel_pw.h"
#include "parallel_global.h"
//#include "../src_algorithms/mymath.h"
//#include "fftw.h"

//#include <unistd.h>

FFT::FFT()
{
	//TITLE("FFT","FFT");
	this->plan_nx = 0;
	this->plan_ny = 0;
	this->plan_nz = 0;
	this->nxx = 0;
	this->scale_xyz = 1;
	this->FFTWsetupwasdone = 0;//mohan modify 2007-11-12
	this->aux4plan = new complex<double>[1];

#ifdef __MPI
	this->plane = new int[1];
	this->aux = new complex<double>[1];
	
	
	this->sentc = new int[1];
	this->sdis = new int[1];
	
	this->recvc = new int[1];
	this->rdis = new int[1];
	
	this->sum = new int[1];
#endif
}

FFT::~FFT()
{
#ifdef __MPI
	delete[] plane;
	delete[] aux;
	delete[] aux4plan;

	
	delete[] sentc;
	delete[] sdis;
	
	delete[] recvc;
	delete[] rdis;
	
	delete[] sum;
#endif
}

void FFT::FFT3D(complex<double> *psi,const int sign)
{
	/*
	cout << "\n\n do fft3d()  "; // << psi[0] << endl;
	cout << "psi size = " << sizeof(psi) << " " << sizeof( psi[0] ) ;
	cout.setf(ios::right);
	cout.setf(ios::scientific);
	cout.precision(15);
	cout << "\n before FFTW, dim: " << plan_nx <<"*"<< plan_ny <<"*"<< plan_nz <<"="<< plan_nx*plan_ny*plan_nz  << "|" << this->nxx << endl;
	for(int i=0; i<3 ; i++) cout<<"\n"<<setw(25)<<psi[i].real()<<setw(25)<<psi[i].imag();
	cout << "\n ... " ;
	for(int i=nxx-3; i<(nxx) ; i++) cout<<"\n"<<setw(25)<<psi[i].real()<<setw(25)<<psi[i].imag();
	cout << endl;
	*/

	timer::tick("FFT","FFT3D");

#ifdef __MPI
	P3DFFT(psi,sign);

#else
	SFFT3D(psi,sign);
#endif

	timer::tick("FFT","FFT3D");

	/*
	cout << "\n\n after FFTW:  \n ";
	for(int i=0; i<3 ; i++) 			cout<<"\n"<<setw(25)<<psi[i].real()<<setw(25)<<psi[i].imag();
	cout << "\n ... " ;
	for(int i=nxx-3; i<(nxx) ; i++)  cout<<"\n"<<setw(25)<<psi[i].real()<<setw(25)<<psi[i].imag();
	cout << "\n ---------------------- :FFT3D sign: " << sign << endl;
	sleep(5);
	*/

	return;
}

/*
void  FFT::FFT3D(double *psi, const int sign)
{
//    cout << "\n do nothing in fft3d() ";
}

void FFT::FFT3D(matrix &psi, const int sign)
{
//    cout << "\n do nothing in fft3d() ";
}
*/



#ifndef __MPI

void FFT::setupFFT3D(const int nx, const int ny, const int nz)
{
	if(test_fft) TITLE("FFT","setupFFT3D");

	this->plan_nx = nx;
	this->plan_ny = ny;
	this->plan_nz = nz;
	this->nxx = nx * ny * nz;
	this->scale_xyz = 1.0 / (double)(this->nxx);  //mohan 2007-6-23
}

void FFT::setupFFT3D_2()
{
	//timer::tick("FFT","setupFFT3D_2");
	
#if defined __FFTW2

	this->plus_plan  = fftw3d_create_plan
	(
		this->plan_nx,this->plan_ny,this->plan_nz,
		FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE | FFTW_USE_WISDOM
	);
	this->minus_plan = fftw3d_create_plan
	(
		this->plan_nx,this->plan_ny,this->plan_nz,
		FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE | FFTW_USE_WISDOM
	);

#elif defined __FFTW3

	delete[] aux4plan; 
	this->aux4plan = new complex<double>[this->nxx];

	//fftw_complex *psiout;
	//psiout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->nxx );
	//std::complex<double>* psiout = new std::complex<double>[this->nxx];
	//fftw_complex *psi2 = reinterpret_cast<fftw_complex*>(psi);
	
	//cout << "\n sign =  1, K->R space \n" ;
	this->plus_plan  = fftw_plan_dft_3d
	(
	//this->plan_nx,this->plan_ny,this->plan_nz, reinterpret_cast<fftw_complex*>(psi),  psiout,
	//this->plan_nx,this->plan_ny,this->plan_nz, reinterpret_cast<fftw_complex*>(psi),  reinterpret_cast<fftw_complex*>(psi) ,
		this->plan_nx,this->plan_ny,this->plan_nz, reinterpret_cast<fftw_complex*>(aux4plan), reinterpret_cast<fftw_complex*>(aux4plan),
		FFTW_BACKWARD, FFTW_ESTIMATE // FFTW_ESTIMATE //FFTW_MEASURE
	);
	//cout << "\n sign = -1, R->K space \n";
	this->minus_plan = fftw_plan_dft_3d
	(
		this->plan_nx,this->plan_ny,this->plan_nz, reinterpret_cast<fftw_complex*>(aux4plan), reinterpret_cast<fftw_complex*>(aux4plan),
		FFTW_FORWARD, FFTW_ESTIMATE  // FFTW_ESTIMATE FFTW_MEASURE
	);

#endif

	if (!this->plus_plan || !this->minus_plan)
	{
		cout << "\nCan't create plans for FFTW in setupFFT3D()\n\n";
	}
	this->FFTWsetupwasdone = true;

	//timer::tick("FFT","setupFFT3D_2");

	return;
}

void FFT::SFFT3D(complex<double> *psi, const int sign)
{
	if(!FFTWsetupwasdone) 
	{
		WARNING("FFT3D","init setupFFT3D_2");
		this->setupFFT3D_2();
	}
	
	if (sign == 1)// K-->R space
	{

#if defined __FFTW2
		fftwnd(this->plus_plan, 1, (FFTW_COMPLEX *)psi, 1, 0, NULL, 0, 0);
#elif defined __FFTW3
		//fftw_execute( this->plus_plan);
		fftw_execute_dft( this->plus_plan, (FFTW_COMPLEX *)psi, (FFTW_COMPLEX *)psi);
#endif

	}
	else if (sign == -1)// R-->K space
	{

#if defined __FFTW2
		fftwnd( this->minus_plan, 1, (FFTW_COMPLEX *)psi, 1, 0, NULL, 0, 0);
#elif defined __FFTW3
		//fftw_execute( this->minus_plan);
		fftw_execute_dft( this->minus_plan, (FFTW_COMPLEX *)psi, (FFTW_COMPLEX *)psi);
#endif

		for (int i = 0; i < this->nxx; i++)
		{
			psi[i] *= this->scale_xyz;
			//psi[i].y *= scale;
		}
	}

	return;
}


#elif defined __MPI

void FFT::setup_MPI_FFT3D(const int nx, const int ny, const int nz, const int nxx_in,const bool in_pool2)
{
	//timer::tick("FFT","Setup_MPI_FFT3D");
	
	if(test_fft) TITLE("FFT","setup_MPI_FFT3D");
	this->plan_nx = nx;
	this->plan_ny = ny;
	this->plan_nz = nz;
	this->nxx = nxx_in;
	this->scale_xy = 1.0 / (double)(nx * ny);
	this->scale_z = 1.0 / (double)(nz);
	this->in_pool = in_pool2;

	if (in_pool)
	{
		this->rank_use = RANK_IN_POOL;
		this->nproc_use = NPROC_IN_POOL;
	}
	else
	{
		this->rank_use = MY_RANK;
		this->nproc_use = NPROC;
	}

	/*
	cout << "\n nx = " << this->plan_nx;
	cout << "\n ny = " << this->plan_ny;
	cout << "\n nz = " << this->plan_nz;
	cout << "\n nxx = " << this->nxx;
	*/

	delete[] plane;
	this->plane = new int[nx];// number of y-z plane
	ZEROS(plane, nx);

	// Searching in all y-z planes
	int i;

	for (int ip=0; ip<nproc_use; ip++)
	{
		for (int is=0; is<nst_per[ip]; is++)
		{
			// find out the stick number.
			int ns = is + this->st_start[ip];

			// find out "x" of FFT grid.
			int occupied = this->ismap[ns]%nx; // mohan 2009-11-09 ny-->nx

			assert(occupied < nx);// mohan 2009-11-09
			this->plane[occupied] = 1;
		}
	}
	/*
	if (!in_pool)
	{
		MPI_Finalize();
		exit(0);
	}
	*/

	delete[] aux;
	this->aux = new complex<double>[this->nxx];

	delete[] sentc;
	delete[] sdis;
	delete[] rdis;
	delete[] recvc;
	delete[] sum;
	
	this->sentc = new int[nproc_use];
	this->sdis = new int[nproc_use];
	this->rdis = new int[nproc_use];
	this->recvc = new int[nproc_use];
	this->sum = new int[nproc_use];

	//---------------------------------------------
	// sum : starting plane of FFT box.
	//---------------------------------------------
	// the first processor:
	this->sum[0] = 0;
	// > first processor
	for (i = 1;i < nproc_use;i++) this->sum[i] = this->sum[i-1] + this->npps[i-1];

	// Each processor has a set of full sticks,
	// 'rank_use' processor send a piece(npps[i]) of these sticks(nst_per[rank_use])
	// to all the other processors in this pool
	for (i = 0;i < nproc_use;i++) this->sentc[i] = this->nst_per[rank_use] * this->npps[i];


	// Each processor in a pool send a piece of each stick(nst_per[i]) to
	// other processors in this pool
	// rank_use processor receive datas in npps[rank_p] planes.
	for (i = 0;i < nproc_use;i++) this->recvc[i] = this->nst_per[i] * this->npps[rank_use];


	// sdis record the starting 'sentc' position in each processor.
	this->sdis[0] = 0;
	for (i = 1;i < nproc_use;i++) this->sdis[i] = this->sdis[i-1] + this->sentc[i-1];


	// rdis record the starting 'recvc' position
	this->rdis[0] = 0;
	for (i = 1;i < nproc_use;i++) this->rdis[i] = this->rdis[i-1] + this->recvc[i-1];


#if defined __FFTW3

	int np = this->npps[rank_use];
	int npy = np * this->plan_ny;
	int ns = this->nst_per[rank_use];
	/*
	cout << " rank_use = " << rank_use << "\n";
	cout << " ns  = " << ns  << "\n";
	cout << " np  = " << np  << "\n";
	cout << " npy = " << npy << "\n";
	*/

	delete[] aux4plan; 
	this->aux4plan = new complex<double>[this->nxx];

	//		                             (rank, *n,			howmany, *in,				  			  *inembed, istride, idist,
	this->planplus_x  = fftw_plan_many_dft(  1, &plan_nx,	npy,	 (fftw_complex *)aux4plan, 		  &plan_nx, npy,     1,
			(fftw_complex *)aux4plan, 	 &plan_nx, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
	//      *out,				 		 *onembed, ostride,	odist,	 sign, 		 	unsigned flags

	//		                             (rank, *n,			howmany, *in,						 	  *inembed, istride, idist,
	this->planplus_y  = fftw_plan_many_dft(  1, &plan_ny,	np,		 (fftw_complex*)&aux4plan[i*npy], &plan_ny, np,      1,
		(fftw_complex*)&aux4plan[i*npy], &plan_ny, np,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
	//	*out, 							 *onembed, ostride,	odist,	 sign, 		 	unsigned flags

	//		                             (rank, *n,			howmany, *in,			  				  *inembed, istride, idist,
	this->planminus_x = fftw_plan_many_dft(  1, &plan_nx, 	npy, 	 (fftw_complex *)aux4plan,		  &plan_nx, npy,     1,
			(fftw_complex *)aux4plan, 	 &plan_nx, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
	//      *out, 				 		 *onembed, ostride,	odist,	 sign, 		 	unsigned flags

	//		                             (rank, *n,			howmany, *in,							  *inembed, istride, idist,
	this->planminus_y = fftw_plan_many_dft(  1, &plan_ny, 	np,		 (fftw_complex*)&aux4plan[i*npy], &plan_ny, np,     1,
		(fftw_complex*)&aux4plan[i*npy], &plan_ny, np,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
	//	*out, 							 *onembed, ostride,	odist,	 sign, 		 	unsigned flags

	//		                             (rank, *n,			howmany, *in,					    	  *inembed, istride, idist,
	this->planplus_z  = fftw_plan_many_dft(  1, &plan_nz,	ns,		 (fftw_complex *)aux4plan,  	  &plan_nz, 1,       plan_nz,
			(fftw_complex *)aux4plan, 	 &plan_nz, 1,		plan_nz, FFTW_BACKWARD,	FFTW_MEASURE   );
	//      *out, 						 *onembed, ostride,  odist,	 sign, 		 	unsigned flags

	//		                             (rank, *n,			howmany, *in,					    	  *inembed, istride, idist,
	this->planminus_z = fftw_plan_many_dft(  1, &plan_nz,	ns,		 (fftw_complex *)aux4plan,  	  &plan_nz, 1,       plan_nz,
			(fftw_complex *)aux4plan, 	 &plan_nz, 1,	    plan_nz, FFTW_FORWARD,	FFTW_MEASURE   );
	//      *out, 						 *onembed, ostride,  odist,	 sign, 		 	unsigned flags

#elif defined __FFTW2

	this->planplus_x = fftw_create_plan(plan_nx, FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE |
	                                  FFTW_THREADSAFE | FFTW_USE_WISDOM);

	this->planplus_y = fftw_create_plan(plan_ny, FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE |
	                                  FFTW_THREADSAFE | FFTW_USE_WISDOM);

	this->planplus_z = fftw_create_plan(plan_nz, FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE |
	                                  FFTW_THREADSAFE | FFTW_USE_WISDOM);

	this->planminus_x = fftw_create_plan(plan_nx, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE |
	                                   FFTW_THREADSAFE | FFTW_USE_WISDOM);

	this->planminus_y = fftw_create_plan(plan_ny, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE |
	                                   FFTW_THREADSAFE | FFTW_USE_WISDOM);

	this->planminus_z = fftw_create_plan(plan_nz, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE |
	                                   FFTW_THREADSAFE | FFTW_USE_WISDOM);

#endif

	if (!planplus_x ||!planplus_y ||!planplus_z || !planminus_x || !planminus_y|| !planminus_z )
	{
		cout << "\nCan't create plans for FFTW in setupFFT3D()\n\n";
		QUIT();
	}

	this->FFTWsetupwasdone = 1;

	//timer::tick("FFT","Setup_MPI_FFT3D");

	if(test_fft)cout << "\n FFTW setup done";
	return;
}

// parallel 3D FFT.
void FFT::P3DFFT(complex<double> *psi, const int sign)
{
	
	//timer::tick("FFT","P3DFFT_Init");

	ZEROS(this->aux, this->nxx);

	// number of z in this cpu.
	const int npps_now = this->npps[rank_use];

	//timer::tick("FFT","P3DFFT_Init");

	// G --> real space
	if (sign == 1)
	{
		this->fftz(psi, sign, this->aux);

		// scatter psi of different cpu.
		// the result is recorded in aux.
		this->scatter(psi, sign);
		
		//timer::tick("FFT","P3DFFT_Map+");
		ZEROS(psi, this->nxx);
		int ii = 0;
		for (int ip = 0;ip < nproc_use;ip++)
		{
			for (int is = 0;is < nst_per[ip];is++)
			{
				// st_start : start stick position
				// in processor i of this pool
				// ns: which (x,y)
				int ir = this->ismap[is + this->st_start[ip]];

				// npps_now: number of z in this cpu.
				for (int k = 0;k < npps_now;k++)
				{
					psi[ir*npps_now+k] = this->aux[ii*npps_now+k];
				}

				ii++;
			}
		}
		//timer::tick("FFT","P3DFFT_Map+");

		this->fftxy(psi, sign);
	}
	else if (sign == -1)
	{
		this->fftxy(psi, sign);
		int sticknow = 0;
		
		//timer::tick("FFT","P3DFFT_Map-");
		for (int ip = 0;ip < nproc_use;ip++)
		{
			for (int j = 0;j < nst_per[ip];j++)
			{
				// use stick number to find (x*n2+y) in FFT grid.
				const int ir = this->ismap[j + st_start[ip]];

				for (int iz=0; iz < npps_now; iz++)
				{
					this->aux[sticknow*npps_now+iz] = psi[ir*npps_now+iz];
				}

				sticknow++;
			}
		}
		//timer::tick("FFT","P3DFFT_Map-");

		this->scatter(psi, sign);
		this->fftz(aux, sign, psi);
	}

	return;
}

void FFT::fftxy(complex<double> *psi, const int sign)
{
	//timer::tick("FFT","fftxy");

	// Number of z in this cpu. 
	int np = this->npps[rank_use];

	// how many points in plane y-z.
	int npy = np * plan_ny;

	if (sign == 1)
	{

		//cout << "np+  = " << np << endl;
		//cout << "npy+ = " << npy << endl;
		//cout << "plan_nx+ = " << plan_nx << endl;

		for (int i=0; i<this->plan_nx; i++)
		{
//			if (this->plane[i] != 0)
//			{
#if defined __FFTW3
				//		                             (rank, *n,				howmany, *in,						 *inembed, istride, idist,
				//this->planplus_y = fftw_plan_many_dft( 1, &plan_ny, 		np,		 (fftw_complex*)&psi[i*npy], &plan_ny, np,      1,
				//		(fftw_complex*)&psi[i*npy], &plan_ny, np,		1,		 FFTW_BACKWARD,	FFTW_ESTIMATE   );
				//       *out, 						*onembed, ostride,	odist,	 sign, 		 	unsigned flags
				//fftw_execute( this->planplus_y );
				//fftw_destroy_plan ( this->planplus_y );
				fftw_execute_dft( this->planplus_y, (fftw_complex*)&psi[i*npy], (fftw_complex*)&psi[i*npy] );
				//fftw_execute_dft( this->planminus_y, (fftw_complex*)&psi[i*npy], (fftw_complex*)&psi[i*npy] );

#elif defined __FFTW2
				fftw(planplus_y, np, (fftw_complex*)&psi[i*npy], np, 1, NULL, 0, 0);
				//fftw(planminus_y, np, (fftw_complex*)&psi[i*npy], np, 1, NULL, 0, 0);
#endif
//			}
		}

#if defined __FFTW3
				//		                            (rank, *n,				howmany, *in,				  *inembed, istride, idist,
				//this->planplus_x  = fftw_plan_many_dft( 1, &plan_nx, 		npy,	 (fftw_complex *)psi, &plan_nx, npy,     1,
				//		(fftw_complex *)psi, &plan_nx, npy,		1,		 FFTW_BACKWARD,	FFTW_ESTIMATE   );
				//      *out,				 *onembed, ostride,	odist,	 sign, 		 	unsigned flags
				//fftw_execute( this->planplus_x );
				//fftw_destroy_plan ( this->planplus_x );
				fftw_execute_dft( this->planplus_x, (fftw_complex *)psi, (fftw_complex *)psi);

#elif defined __FFTW2
				fftw(planplus_x, npy, (fftw_complex *)psi, npy, 1, NULL, 0, 0);
#endif

	}
	else if (sign == -1)
	{

		//cout << "np-  = " << np << endl;
		//cout << "npy- = " << npy << endl;
		//cout << "plan_nx- = " << plan_nx << endl;

#if defined __FFTW3
				//		                               (rank, *n,				howmany, *in,			  *inembed, istride, idist,
				//this->planminus_x  = fftw_plan_many_dft( 1, &plan_nx, 		npy, (fftw_complex *)psi, &plan_nx, npy,     1,
				//		(fftw_complex *)psi, &plan_nx, npy,		1,		 FFTW_FORWARD,	FFTW_ESTIMATE   );
				//      *out, 				 *onembed, ostride,	odist,	 sign, 		 	unsigned flags
				//fftw_execute( this->planminus_x );
				//fftw_destroy_plan ( this->planminus_x );
				fftw_execute_dft( this->planminus_x, (fftw_complex *)psi, (fftw_complex *)psi);

#elif defined __FFTW2
				fftw(planminus_x, npy, (fftw_complex *)psi, npy, 1, NULL, 0, 0);
#endif

		for (int i = 0;i <this->plan_nx;i++)
		{
//			if (this->plane[i] != 0)
//			{

#if defined __FFTW3
				//		                              (rank, *n,	   howmany, *in,						*inembed, istride, idist,
				//this->planminus_y = fftw_plan_many_dft( 1, &plan_ny, 		np,	(fftw_complex*)&psi[i*npy], &plan_ny, np,     1,
				//		(fftw_complex*)&psi[i*npy], &plan_ny, np,		1,		 FFTW_FORWARD,	FFTW_ESTIMATE   );
				//       *out, 						*onembed, ostride,	odist,	 sign, 		 	unsigned flags
				//fftw_execute( this->planminus_y );
				//fftw_destroy_plan ( this->planminus_y );
				fftw_execute_dft( this->planminus_y, (fftw_complex*)&psi[i*npy], (fftw_complex*)&psi[i*npy]);

#elif defined __FFTW2
				fftw(planminus_y, np, (fftw_complex *) &psi[i*npy], np, 1, NULL, 0, 0);
//				fftw(planminus_y, np, (fftw_complex *) &psi[i*npy], np, plan_nx, NULL, 0, 0);
#endif
//			}
		}

		for (int i = 0;i < plan_nx*npy;i++)
		{
			psi[i] *= scale_xy;
		}
	}
	
	//timer::tick("FFT","fftxy");
	return;
}



void FFT::fftz(complex<double> *psi_in, const int sign, complex<double> *psi_out)
{
	//timer::tick("FFT","fftz");
	// number of sticks in this process.
#if defined __FFTW2
	int ns = this->nst_per[rank_use];
#endif

	if (sign == 1)
	{
		// cout << "ns+ = " << ns << endl;
		// only do ns * plan_nz(element number) fft .

#if defined __FFTW3
		//		                             (rank, *n,				 howmany, *in,					  *inembed, istride, idist,
		//this->planplus_z = fftw_plan_many_dft( 1, &plan_nz, 		 ns,	  (fftw_complex *)psi_in, &plan_nz, 1,       plan_nz,
		//		(fftw_complex *)psi_in, &plan_nz, 1,		 plan_nz, FFTW_BACKWARD,	FFTW_ESTIMATE   );
		//      *out, 					*onembed, ostride, odist,	  sign, 		 	unsigned flags
		//fftw_execute( this->planplus_z );
		//fftw_destroy_plan ( this->planplus_z );
		fftw_execute_dft( this->planplus_z, (fftw_complex *)psi_in, (fftw_complex *)psi_in );

#elif defined __FFTW2
		//this->planplus_z = fftw_create_plan(plan_nz, FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE |
		//                                  FFTW_THREADSAFE | FFTW_USE_WISDOM);
		fftw(planplus_z, ns, (fftw_complex *)psi_in, 1, plan_nz, NULL, 0, 0);
#endif

		for (int i = 0;i < this->nxx;i++)
		{
			psi_out[i] = psi_in[i];
		}
	}
	else if (sign == -1)
	{

		// cout << "ns- = " << ns << endl;

#if defined __FFTW3
		//		                              (rank, *n,			 howmany, *in,					  *inembed, istride, idist,
		//this->planminus_z = fftw_plan_many_dft( 1, &plan_nz, 		 ns,	  (fftw_complex *)psi_in, &plan_nz, 1,       plan_nz,
		//		(fftw_complex *)psi_in, &plan_nz, 1,	   plan_nz, FFTW_FORWARD,	FFTW_ESTIMATE   );
		//       *out, 					*onembed, ostride, odist,	  sign, 		 	unsigned flags
		//fftw_execute( this->planminus_z );
		//fftw_destroy_plan ( this->planminus_z );
		fftw_execute_dft( this->planminus_z, (fftw_complex *)psi_in, (fftw_complex *)psi_in );

#elif defined __FFTW2
		fftw(planminus_z, ns, (fftw_complex *)psi_in, 1, plan_nz, NULL, 0, 0);
#endif

		for (int i = 0;i < this->nxx;i++)
		{
			psi_out[i] = psi_in[i] * scale_z;
		}
	}

	//timer::tick("FFT","fftz");
	return;
}



void FFT::scatter(complex<double> *psi, int sign)
{
	//timer::tick("FFT","scatter");

	int ns = nst_per[rank_use];
	int nz = this->plan_nz;

	if (nproc_use == 1)
	{
		return;
	}

	if (sign == 1)
	{
		for (int i = 0;i < nproc_use;i++)
		{
			int npps_now = this->npps[i];

			for (int j = 0;j < ns;j++)
			{
				for (int k = 0;k < npps_now;k++)
				{
					//------------------------------------------------------
					// i : distinguish different cpu in this pool
					// j : distinguish different stick in this pool
					// k : each plane in cpu i.
					// j*nz : starting position of stick j.
					// sum[i] :  starting position of x-y planes,
					// 			 these planes(npps_now) belongs to cpu i.
					// npps_now :  number of x-y planes in cpu i.
					// sdis: starting position of 'sendc'
					//
					// prepare for communication between process:
					// cutting sticks into pieces , according to
					// different plane number in different cpu.
					//------------------------------------------------------
					psi[sdis[i] + j*npps_now + k] = this->aux[ j*nz + sum[i] + k];
				}
			}
		}

		ZEROS(this->aux, nxx);

		if (in_pool)
		{
			// send buffer ==> psi
			// receive buffer ==> aux 
			MPI_Alltoallv(psi, sentc, sdis, mpicomplex, this->aux, recvc, rdis, mpicomplex, POOL_WORLD);
		}
		else
		{
			MPI_Alltoallv(psi, sentc, sdis, mpicomplex, this->aux, recvc, rdis, mpicomplex, MPI_COMM_WORLD);
		}
	}
	else if (sign == -1)
	{
		if (in_pool)
		{
			MPI_Alltoallv(this->aux, recvc, rdis, mpicomplex, psi, sentc, sdis, mpicomplex, POOL_WORLD);
		}
		else
		{
			MPI_Alltoallv(this->aux, recvc, rdis, mpicomplex, psi, sentc, sdis, mpicomplex, MPI_COMM_WORLD);
		}

		ZEROS(this->aux, nxx);

		for (int i = 0;i < nproc_use;i++)
		{
			int nplane = npps[i];

			for (int j = 0;j < ns;j++)
			{
				for (int k = 0;k < nplane;k++)
				{
					this->aux[j*nz+sum[i] + k] = psi[sdis[i] + j * nplane + k];
				}
			}
		}
	}

	//timer::tick("FFT","scatter");
	return;
}
#endif // __MPI
