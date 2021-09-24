#include "global.h"
#include "hamilt.h"
#include "diago_cg.h"
#include "diago_david.h"
#include "diago_cg_gpu.h"
#include "cufft.h"



Hamilt::Hamilt() {}
Hamilt::~Hamilt() {}

// in tools.h

__global__ void hamilt_cast_d2f(float *dst, double *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i] = __double2float_rn(src[i]);
    }
}

__global__ void hamilt_cast_f2d(double *dst, float *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i] = (double)(src[i]);
    }
}

__global__ void hamilt_cast_d2f(float2 *dst, double2 *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i].x = __double2float_rn(src[i].x);
        dst[i].y = __double2float_rn(src[i].y);
    }
}

__global__ void hamilt_cast_f2d(double2 *dst, float2 *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i].x = (double)(src[i].x);
        dst[i].y = (double)(src[i].y);
    }
}


void Hamilt::diagH_pw(
    const int &istep,
    const int &iter,
    const int &ik,
    const double *precondition,
    double &avg_iter)
{
	ModuleBase::TITLE("Hamilt","diagH_pw");
    ModuleBase::timer::tick("Hamilt", "diagH_pw");
    double avg = 0.0;

	// set ik0 because of mem_saver.
	// if mem_saver is not used, ik0=ik, as usual.
	// but if mem_saver is used, ik0=0.
	int ik0 = ik;

	if(GlobalV::CALCULATION=="nscf" && GlobalC::wf.mem_saver==1)
	{
		if(GlobalV::BASIS_TYPE=="pw")
		{
			// generate PAOs first, then diagonalize to get
			// inital wavefunctions.
			GlobalC::wf.diago_PAO_in_pw_k2(ik, GlobalC::wf.evc[0]);
		}
#ifdef __LCAO
		else if(GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			GlobalC::wf.LCAO_in_pw_k(ik, GlobalC::wf.wanf2[0]);
		}
#endif
		ik0 = 0;
	}

    if(GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
		if(GlobalV::KS_SOLVER=="lapack")
		{
			assert(GlobalV::NLOCAL >= GlobalV::NBANDS);
        	this->diagH_subspace(
				ik,
				GlobalV::NLOCAL,
				GlobalV::NBANDS,
				GlobalC::wf.wanf2[ik0],
				GlobalC::wf.evc[ik0],
				GlobalC::wf.ekb[ik]);
		}
		else
		{
			GlobalV::ofs_warning << " The diago_type " << GlobalV::KS_SOLVER
				<< " not implemented yet." << std::endl; //xiaohui add 2013-09-02
            ModuleBase::WARNING_QUIT("Hamilt::diago","no implemt yet.");
		}
    }
    else
    {
        int ntry = 0;
        int notconv = 0;
        do
        {
	   		if(GlobalV::KS_SOLVER=="cg")
            {
				// qian change it, because it has been executed in diago_PAO_in_pw_k2
                if ( iter > 1 || istep > 1 ||  ntry > 0)
                {
                    this->diagH_subspace(
						ik,
						GlobalV::NBANDS,
						GlobalV::NBANDS,
						GlobalC::wf.evc[ik0],
						GlobalC::wf.evc[ik0],
						GlobalC::wf.ekb[ik]);

                    avg_iter += 1.0;
                }

                // Diago_CG_CUDA<float, float2> f_cg_gpu;
                Diago_CG_CUDA<double, double2> d_cg_gpu;
                
				bool reorder = true;

				double2 *d_wf_evc;
                double *d_wf_ekb;
                int DIM_CG_GPU = GlobalC::kv.ngk[ik];
                int DIM_CG_GPU2 = GlobalC::wf.npwx * GlobalV::NPOL;
                double *d_precondition;

                float2 *f_wf_evc;
                float *f_wf_ekb;
                // int DIM_CG_GPU = GlobalC::kv.ngk[ik];
                // int DIM_CG_GPU2 = GlobalC::wf.npwx * GlobalV::NPOL;
                float *f_precondition;

				if(GlobalV::NPOL==1)
				{
                    CHECK_CUDA(cudaMalloc((void**)&d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(double2)));
                    CHECK_CUDA(cudaMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double)));
                    CHECK_CUDA(cudaMalloc((void**)&d_precondition, DIM_CG_GPU * sizeof(double)));

                    CHECK_CUDA(cudaMemcpy(d_wf_evc, GlobalC::wf.evc[ik0].c, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(double2), cudaMemcpyHostToDevice));
                    // CHECK_CUDA(cudaMemcpy(d_wf_ekb, wf.ekb[ik], NBANDS * sizeof(double), cudaMemcpyHostToDevice));
                    CHECK_CUDA(cudaMemcpy(d_precondition, precondition, DIM_CG_GPU * sizeof(double), cudaMemcpyHostToDevice));
                    
                    // cast to float
                    // if(istep < 3)
                    // {
                    //     CHECK_CUDA(cudaMalloc((void**)&f_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(float2)));
                    //     CHECK_CUDA(cudaMalloc((void**)&f_wf_ekb, GlobalV::NBANDS * sizeof(float)));
                    //     CHECK_CUDA(cudaMalloc((void**)&f_precondition, DIM_CG_GPU * sizeof(float)));
                    //     int thread = 512;
                    //     int block = GlobalV::NBANDS * GlobalC::wf.npwx / thread + 1;
                    //     int block2 = GlobalV::NBANDS / thread + 1;
                    //     int block3 = DIM_CG_GPU / thread + 1;
                    //     hamilt_cast_d2f<<<block, thread>>>(f_wf_evc, d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx);
                    //     hamilt_cast_d2f<<<block3, thread>>>(f_precondition, d_precondition, DIM_CG_GPU);
                    //     cufftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_C2C);
                    //     f_cg_gpu.diag(f_wf_evc, f_wf_ekb, DIM_CG_GPU, GlobalC::wf.npwx,
                    //         GlobalV::NBANDS, f_precondition, GlobalV::ETHR,
                    //         GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);
                    //     hamilt_cast_f2d<<<block, thread>>>(d_wf_evc, f_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx);
                    //     hamilt_cast_f2d<<<block2, thread>>>(d_wf_ekb, f_wf_ekb, GlobalV::NBANDS);
                    // }
                    // else
                    // {
                    cufftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
                    d_cg_gpu.diag(d_wf_evc, d_wf_ekb, DIM_CG_GPU, GlobalC::wf.npwx,
                        GlobalV::NBANDS, d_precondition, GlobalV::ETHR,
                        GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);
                    // }
                    
                    cufftDestroy(GlobalC::UFFT.fft_handle);
                    // to cpu
                    CHECK_CUDA(cudaMemcpy(GlobalC::wf.evc[ik0].c, d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(double2), cudaMemcpyDeviceToHost));
                    CHECK_CUDA(cudaMemcpy(GlobalC::wf.ekb[ik], d_wf_ekb, GlobalV::NBANDS * sizeof(double), cudaMemcpyDeviceToHost));

                    CHECK_CUDA(cudaFree(d_wf_evc));
                    CHECK_CUDA(cudaFree(d_wf_ekb));
                    CHECK_CUDA(cudaFree(d_precondition));
				}
                // comment this to debug.
                
				else
				{
                    // to gpu
                    CHECK_CUDA(cudaMalloc((void**)&d_wf_evc, GlobalV::NBANDS * DIM_CG_GPU2 * sizeof(double2)));
                    CHECK_CUDA(cudaMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double)));
                    CHECK_CUDA(cudaMalloc((void**)&d_precondition, DIM_CG_GPU2 * sizeof(double)));

                    CHECK_CUDA(cudaMemcpy(d_wf_evc, GlobalC::wf.evc[ik0].c, GlobalV::NBANDS * DIM_CG_GPU2 * sizeof(double2), cudaMemcpyHostToDevice));
                    // CHECK_CUDA(cudaMemcpy(d_wf_ekb, GlobalC::wf.ekb[ik], GlobalV::NBANDS * sizeof(double), cudaMemcpyHostToDevice));
                    CHECK_CUDA(cudaMemcpy(d_precondition, precondition, DIM_CG_GPU2 * sizeof(double), cudaMemcpyHostToDevice));
                    // do things
                    cufftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);

                    d_cg_gpu.diag(d_wf_evc, d_wf_ekb, DIM_CG_GPU2, DIM_CG_GPU2,
                        GlobalV::NBANDS, d_precondition, GlobalV::ETHR,
                        GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);
                    cufftDestroy(GlobalC::UFFT.fft_handle);

                    // to cpu
                    CHECK_CUDA(cudaMemcpy(GlobalC::wf.evc[ik0].c, d_wf_evc, GlobalV::NBANDS * DIM_CG_GPU2 * sizeof(double2), cudaMemcpyDeviceToHost));
                    CHECK_CUDA(cudaMemcpy(GlobalC::wf.ekb[ik], d_wf_ekb, GlobalV::NBANDS * sizeof(double), cudaMemcpyDeviceToHost));

                    CHECK_CUDA(cudaFree(d_wf_evc));
                    CHECK_CUDA(cudaFree(d_wf_ekb));
                    CHECK_CUDA(cudaFree(d_precondition));
				}
                
				// P.S. : nscf is the flag about reorder.
				// if diagH_subspace is done once,
				// we don't need to reorder the eigenvectors order.
				// if diagH_subspace has not been called,
				// we need to reorder the eigenvectors.
            }
	   		else if(GlobalV::KS_SOLVER=="dav")
        	{
				Diago_David david;
				if(GlobalV::NPOL==1)
				{
					david.diag(GlobalC::wf.evc[ik0], GlobalC::wf.ekb[ik], GlobalC::kv.ngk[ik],
						GlobalV::NBANDS, precondition, GlobalV::DIAGO_DAVID_NDIM,
				 		GlobalV::ETHR, GlobalV::DIAGO_CG_MAXITER, notconv, avg);
				}
				else
				{
					david.diag(GlobalC::wf.evc[ik0], GlobalC::wf.ekb[ik], GlobalC::wf.npwx*GlobalV::NPOL,
						GlobalV::NBANDS, precondition, GlobalV::DIAGO_DAVID_NDIM,
						GlobalV::ETHR, GlobalV::DIAGO_CG_MAXITER, notconv, avg);
				}
        	}
        	else
        	{
				ModuleBase::WARNING_QUIT("calculate_bands","Check GlobalV::KS_SOLVER !");
        	}
            avg_iter += avg;
            ++ntry;
        }
        while ( this->test_exit_cond(ntry, notconv) );

        if ( notconv > max(5, GlobalV::NBANDS/4) )
        {
            std::cout << "\n notconv = " << notconv;
            std::cout << "\n Hamilt::diago', too many bands are not converged! \n";
        }
    }

	ModuleBase::timer::tick("Hamilt","diagH_pw");
    return;
}


bool Hamilt::test_exit_cond(const int &ntry, const int &notconv)
{
    //================================================================
    // If this logical function is true, need to do diagH_subspace
	// and cg again.
    //================================================================

	bool scf = true;
	if(GlobalV::CALCULATION=="nscf") scf=false;

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ( (!scf && (notconv > 0)) );

    // if self consistent calculation, if not converged > 5,
    // using diagH_subspace and cg method again. ntry++
    const bool f3 = ( ( scf && (notconv > 5)) );
    return  ( f1 && ( f2 || f3 ) );
}

void Hamilt::diagH_subspace(
    const int ik,
    const int nstart,
    const int n_band,
    const ModuleBase::ComplexMatrix &psi,
    ModuleBase::ComplexMatrix &evc,
    double *en)
{
	if(nstart < n_band)
	{
		ModuleBase::WARNING_QUIT("diagH_subspace","nstart < n_band!");
	}

    if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
        this->hpw.diagH_subspace(ik, nstart, n_band, psi, evc, en);
    }
    else
    {
		ModuleBase::WARNING_QUIT("diagH_subspace","Check parameters: GlobalV::BASIS_TYPE. ");
    }
    return;
}



//====================================================================
// calculates eigenvalues and eigenvectors of the generalized problem
// Hv=eSv, with H hermitean matrix, S overlap matrix .
// On output both matrix are unchanged
// LAPACK version - uses both ZHEGV and ZHEGVX
//=====================================================================
void Hamilt::diagH_LAPACK(
	const int nstart,
	const int nbands,
	const ModuleBase::ComplexMatrix &hc,
	const ModuleBase::ComplexMatrix &sc,
	const int ldh, // nstart
	double *e,
	ModuleBase::ComplexMatrix &hvec)
{
    ModuleBase::TITLE("Hamilt","diagH_LAPACK");
	ModuleBase::timer::tick("Hamilt","diagH_LAPACK");

    int lwork=0;
    //========================================
    // int ILAENV();
    // ILAENV returns optimal block size "nb"
    //========================================

    ModuleBase::ComplexMatrix sdum(nstart, ldh);
    ModuleBase::ComplexMatrix hdum;

    sdum = sc;

    const bool all_eigenvalues = (nstart == nbands);

    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);
//  int nb = ILAENV(1,  "ZHETRD", "U", n, -1, -1, -1);
//  int nb = 32;

    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }

	if (nb == 1 || nb >= nstart)
    {
        lwork = 2 * nstart; // mohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }

    std::complex<double> *work = new std::complex<double>[lwork];
	ModuleBase::GlobalFunc::ZEROS(work, lwork);

    //=====================================================================
    // input s and (see below) h are copied so that they are not destroyed
    //=====================================================================

    int info = 0;
    int rwork_dim;
    if (all_eigenvalues)
    {
        rwork_dim = 3*nstart-2;
    }
    else
    {
        rwork_dim = 7*nstart;
    }

    double *rwork = new double[rwork_dim];
    ModuleBase::GlobalFunc::ZEROS( rwork, rwork_dim );

    if (all_eigenvalues)
    {
        //===========================
        // calculate all eigenvalues
        //===========================
        hvec = hc;
        LapackConnector::zhegv(1, 'V', 'U', nstart, hvec , ldh, sdum, ldh, e, work , lwork , rwork, info);
    }
    else
    {
        //=====================================
        // calculate only m lowest eigenvalues
        //=====================================
        int *iwork = new int [5*nstart];
        int *ifail = new int[nstart];

        ModuleBase::GlobalFunc::ZEROS(rwork,7*nstart);
        ModuleBase::GlobalFunc::ZEROS(iwork,5*nstart);
        ModuleBase::GlobalFunc::ZEROS(ifail,nstart);

        hdum.create(nstart, ldh);
        hdum = hc;

    	//=============================
    	// Number of calculated bands
    	//=============================
    	int mm = nbands;

        LapackConnector::zhegvx
        (
                1,      //INTEGER
                'V',    //CHARACTER*1
                'I',    //CHARACTER*1
                'U',    //CHARACTER*1
                nstart, //INTEGER
                hdum,   //COMPLEX*16 array
                ldh,    //INTEGER
                sdum,   //COMPLEX*16 array
                ldh,    //INTEGER
           		0.0,    //DOUBLE PRECISION
                0.0,    //DOUBLE PRECISION
                1,      //INTEGER
                nbands, //INTEGER
                0.0,    //DOUBLE PRECISION
                mm,     //INTEGER
                e,      //DOUBLE PRECISION array
                hvec,   //COMPLEX*16 array
                ldh,    //INTEGER
                work,   //DOUBLE array, dimension (MAX(1,LWORK))
                lwork,  //INTEGER
                rwork , //DOUBLE PRECISION array, dimension (7*N)
                iwork,  //INTEGER array, dimension (5*N)
                ifail,  //INTEGER array, dimension (N)
                info    //INTEGER
        );

        delete[] iwork;
        delete[] ifail;
    }
    delete[] rwork;
    delete[] work;

	ModuleBase::timer::tick("Hamilt","diagH_LAPACK");
    return;
}
