#include "hip/hip_runtime.h"
#include "global.h"
#include "hamilt.h"
// #include "diago_cg.h"
#include "diago_david.h"
#include "diago_cg_hip.h"
#include "hipfft.h"

Hamilt::Hamilt() 
{
}
Hamilt::~Hamilt() 
{
}

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
                
                hipblasDoubleComplex *d_wf_evc;
                double *d_wf_ekb;
                hipblasDoubleComplex *d_vkb_c;

                int nkb = GlobalC::ppcell.nkb;
                if(iter < 0)
                {
                    CHECK_CUFFT(hipfftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_C2C));
                }
                else
                {
                    CHECK_CUFFT(hipfftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
                }

                CHECK_CUDA(hipMalloc((void**)&d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * GlobalV::NPOL * sizeof(hipblasDoubleComplex)));
                CHECK_CUDA(hipMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double)));

                CHECK_CUDA(hipMemcpy(d_wf_evc, 
                            GlobalC::wf.evc[ik0].c, 
                            GlobalV::NBANDS * GlobalC::wf.npwx * GlobalV::NPOL * sizeof(hipblasDoubleComplex), hipMemcpyHostToDevice));

                CHECK_CUDA(hipMalloc((void**)&d_vkb_c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex)));
                CHECK_CUDA(hipMemcpy(d_vkb_c, GlobalC::ppcell.vkb.c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex), hipMemcpyHostToDevice));
                
                if ( iter > 1 || istep > 1 ||  ntry > 0)
                {
                    this->diagH_subspace_cuda(
						ik,
						GlobalV::NBANDS,
						GlobalV::NBANDS,
						d_wf_evc,
						d_wf_evc,
						d_wf_ekb,
                        d_vkb_c);

                    avg_iter += 1.0;
                }

                Diago_CG_CUDA<float, hipblasComplex, float2> f_cg_cuda;
                Diago_CG_CUDA<double, hipblasDoubleComplex, double2> d_cg_cuda;
                
				bool reorder = true;

				
                int DIM_CG_CUDA = GlobalC::kv.ngk[ik];
                int DIM_CG_CUDA2 = GlobalC::wf.npwx * GlobalV::NPOL;
                double *d_precondition;

                hipblasComplex *f_wf_evc;
                float *f_wf_ekb;
                float *f_precondition;
                hipblasComplex *f_vkb_c;

				if(GlobalV::NPOL==1)
				{
                    // CHECK_CUDA(hipMalloc((void**)&d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(hipblasDoubleComplex)));
                    // CHECK_CUDA(hipMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double)));
                    CHECK_CUDA(hipMalloc((void**)&d_precondition, DIM_CG_CUDA * sizeof(double)));
                    
                    // Add d_vkb_c
                    // CHECK_CUDA(hipMalloc((void**)&d_vkb_c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex)));
                    // CHECK_CUDA(hipMemcpy(d_vkb_c, GlobalC::ppcell.vkb.c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex), hipMemcpyHostToDevice));

                    // CHECK_CUDA(hipMemcpy(d_wf_evc, GlobalC::wf.evc[ik0].c, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(hipblasDoubleComplex), hipMemcpyHostToDevice));
                    // CHECK_CUDA(hipMemcpy(d_wf_ekb, wf.ekb[ik], NBANDS * sizeof(double), hipMemcpyHostToDevice));
                    CHECK_CUDA(hipMemcpy(d_precondition, precondition, DIM_CG_CUDA * sizeof(double), hipMemcpyHostToDevice));
                    // cout<<"ITER: "<<iter<<endl;
                    // cout<<"ETHR_now: "<<GlobalV::ETHR<<endl;
                    // cast to float
                    // if(iter < 100)
                    // if(iter == 1 || GlobalV::ETHR > 5e-4)
                    if(iter < 0)
                    {
                        CHECK_CUDA(hipMalloc((void**)&f_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(hipblasComplex)));
                        CHECK_CUDA(hipMalloc((void**)&f_wf_ekb, GlobalV::NBANDS * sizeof(float)));
                        CHECK_CUDA(hipMalloc((void**)&f_precondition, DIM_CG_CUDA * sizeof(float)));

                        // add vkb_c parameter
                        CHECK_CUDA(hipMalloc((void**)&f_vkb_c, GlobalC::wf.npwx*nkb*sizeof(hipblasComplex)));

                        int thread = 512;
                        int block = GlobalV::NBANDS * GlobalC::wf.npwx / thread + 1;
                        int block2 = GlobalV::NBANDS / thread + 1;
                        int block3 = DIM_CG_CUDA / thread + 1;
                        int block4 = GlobalC::wf.npwx*nkb / thread + 1;

                        hipLaunchKernelGGL(hamilt_cast_d2f, dim3(block), dim3(thread), 0, 0, 
                            reinterpret_cast<float2*>(f_wf_evc), reinterpret_cast<double2*>(d_wf_evc), GlobalV::NBANDS * GlobalC::wf.npwx);
                        hipLaunchKernelGGL(hamilt_cast_d2f, dim3(block3), dim3(thread), 0, 0, 
                            reinterpret_cast<float2*>(f_precondition), reinterpret_cast<double2*>(d_precondition), DIM_CG_CUDA);
                        // add vkb_c parameter
                        hipLaunchKernelGGL(hamilt_cast_d2f, dim3(block4), dim3(thread), 0, 0, 
                            reinterpret_cast<float2*>(f_vkb_c), reinterpret_cast<double2*>(d_vkb_c), GlobalC::wf.npwx*nkb);
                        // CHECK_CUFFT(hipfftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_C2C));
                        // cout<<"Do float CG ..."<<endl;
                        f_cg_cuda.diag(f_wf_evc, f_wf_ekb, f_vkb_c, DIM_CG_CUDA, GlobalC::wf.npwx,
                            GlobalV::NBANDS, f_precondition, GlobalV::ETHR,
                            GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);
                        hipLaunchKernelGGL(hamilt_cast_f2d, dim3(block), dim3(thread), 0, 0, 
                            reinterpret_cast<double2*>(d_wf_evc), reinterpret_cast<float2*>(f_wf_evc), GlobalV::NBANDS * GlobalC::wf.npwx);
                        hipLaunchKernelGGL(hamilt_cast_f2d, dim3(block2), dim3(thread), 0, 0, 
                            reinterpret_cast<double2*>(d_wf_ekb), reinterpret_cast<float2*>(f_wf_ekb), GlobalV::NBANDS);

                        CHECK_CUDA(hipFree(f_vkb_c));
                    }
                    else
                    {
                        // cout<<"begin cg!!"<<endl;
                        // CHECK_CUFFT(hipfftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
                        // cout<<"Do double CG ..."<<endl;
                        d_cg_cuda.diag(d_wf_evc, d_wf_ekb, d_vkb_c, DIM_CG_CUDA, GlobalC::wf.npwx,
                            GlobalV::NBANDS, d_precondition, GlobalV::ETHR,
                            GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);
                    }
                    // TODO destroy handle
                    CHECK_CUFFT(hipfftDestroy(GlobalC::UFFT.fft_handle));
                    // to cpu
                    CHECK_CUDA(hipMemcpy(GlobalC::wf.evc[ik0].c, d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(hipblasDoubleComplex), hipMemcpyDeviceToHost));
                    CHECK_CUDA(hipMemcpy(GlobalC::wf.ekb[ik], d_wf_ekb, GlobalV::NBANDS * sizeof(double), hipMemcpyDeviceToHost));

                    // CHECK_CUDA(hipFree(d_wf_evc));
                    // CHECK_CUDA(hipFree(d_wf_ekb));
                    // CHECK_CUDA(hipFree(d_vkb_c))
                    CHECK_CUDA(hipFree(d_precondition));
				}
                // comment this to debug.
                
				else
				{
                    // to gpu
                    // CHECK_CUDA(hipMalloc((void**)&d_wf_evc, GlobalV::NBANDS * DIM_CG_CUDA2 * sizeof(hipblasDoubleComplex)));
                    // CHECK_CUDA(hipMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double)));
                    CHECK_CUDA(hipMalloc((void**)&d_precondition, DIM_CG_CUDA2 * sizeof(double)));

                    // CHECK_CUDA(hipMemcpy(d_wf_evc, GlobalC::wf.evc[ik0].c, GlobalV::NBANDS * DIM_CG_CUDA2 * sizeof(hipblasDoubleComplex), hipMemcpyHostToDevice));
                    // CHECK_CUDA(hipMemcpy(d_wf_ekb, GlobalC::wf.ekb[ik], GlobalV::NBANDS * sizeof(double), hipMemcpyHostToDevice));
                    CHECK_CUDA(hipMemcpy(d_precondition, precondition, DIM_CG_CUDA2 * sizeof(double), hipMemcpyHostToDevice));
                    // do things
                    CHECK_CUFFT(hipfftPlan3d(&GlobalC::UFFT.fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));

                    d_cg_cuda.diag(d_wf_evc, d_wf_ekb, d_vkb_c, DIM_CG_CUDA2, DIM_CG_CUDA2,
                        GlobalV::NBANDS, d_precondition, GlobalV::ETHR,
                        GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);
                    CHECK_CUFFT(hipfftDestroy(GlobalC::UFFT.fft_handle));

                    // to cpu
                    CHECK_CUDA(hipMemcpy(GlobalC::wf.evc[ik0].c, d_wf_evc, GlobalV::NBANDS * DIM_CG_CUDA2 * sizeof(hipblasDoubleComplex), hipMemcpyDeviceToHost));
                    CHECK_CUDA(hipMemcpy(GlobalC::wf.ekb[ik], d_wf_ekb, GlobalV::NBANDS * sizeof(double), hipMemcpyDeviceToHost));

                    // CHECK_CUDA(hipFree(d_wf_evc));
                    // CHECK_CUDA(hipFree(d_wf_ekb));
                    CHECK_CUDA(hipFree(d_precondition));
				}

                CHECK_CUDA(hipFree(d_wf_evc));
                CHECK_CUDA(hipFree(d_wf_ekb));
                CHECK_CUDA(hipFree(d_vkb_c))
                
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

void Hamilt::diagH_subspace_cuda(
    const int ik,
    const int nstart,
    const int n_band,
    const hipblasDoubleComplex* psi,
    hipblasDoubleComplex* evc,
    double *en,
    hipblasDoubleComplex *d_vkb_c)
{
	if(nstart < n_band)
	{
		ModuleBase::WARNING_QUIT("diagH_subspace_cuda","nstart < n_band!");
	}

    if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
        this->hpw.diagH_subspace_cuda(ik, nstart, n_band, psi, evc, en, d_vkb_c);
    }
    else
    {
		ModuleBase::WARNING_QUIT("diagH_subspace_cuda","Check parameters: GlobalV::BASIS_TYPE. ");
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
	const ModuleBase::ComplexMatrix &hc, // nstart * nstart
	const ModuleBase::ComplexMatrix &sc, // nstart * nstart
	const int ldh, // nstart
	double *e,
	ModuleBase::ComplexMatrix &hvec)  // nstart * n_band
{
    ModuleBase::TITLE("Hamilt","diagH_LAPACK");
	ModuleBase::timer::tick("Hamilt","diagH_LAPACK");

    // Print info ... 
    // cout<<"in diagH_lapack"<<endl;
    // cout<<nstart<<endl;
    // cout<<nbands<<endl;
    // cout<<"hc: "<<hc.nr<<" "<<hc.nc<<endl;
    // cout<<"sc: "<<sc.nr<<" "<<sc.nc<<endl;
    // cout<<ldh<<endl;
    // cout<<"hvec: "<<hvec.nr<<" "<<hvec.nc<<endl;

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

/*
void Hamilt::diagH_CUSOLVER(
	const int nstart,
	const int nbands,
	hipblasDoubleComplex* hc,  // nstart * nstart
	hipblasDoubleComplex* sc,  // nstart * nstart
	const int ldh, // nstart
	double *e,
	hipblasDoubleComplex* hvec)  // nstart * n_band
{
    ModuleBase::TITLE("Hamilt","diagH_CUSOLVER");
	ModuleBase::timer::tick("Hamilt","diagH_CUSOLVER");

    // cout<<"diagh rocsolver!!"<<endl;

    // cout<<nstart<<" "<<nbands<<endl;

    hipblasDoubleComplex *sdum;
    hipblasDoubleComplex *hdum;
    CHECK_CUDA(hipMalloc((void**)&sdum, nstart*ldh*sizeof(hipblasDoubleComplex)));
    // CHECK_CUDA(hipMalloc((void**)&hdum, nstart*ldh*sizeof(hipblasDoubleComplex)));

    // sdum = sc;
    CHECK_CUDA(hipMemcpy(sdum, sc, nstart*ldh*sizeof(hipblasDoubleComplex), hipMemcpyDeviceToDevice));
    hipDeviceSynchronize();

    int *device_info;
    int h_meig = 0;
    CHECK_CUDA(hipMalloc((void**)&device_info, sizeof(int)));
    const bool all_eigenvalues = (nstart == nbands);

    // rocsolverDnHandle_t rocsolver_handle;
    // CHECK_CUSOLVER(rocsolverDnCreate(&rocsolver_handle));

    if (all_eigenvalues) // nstart = nbands
    {
        //===========================
        // calculate all eigenvalues
        //===========================
        // hvec = hc;
        // repalce "=" with memcpy!
        CHECK_CUDA(hipMemcpy(hvec, hc, nstart*nbands*sizeof(hipblasDoubleComplex), hipMemcpyDeviceToDevice));
        hipDeviceSynchronize();
        // use rocsolver
        int rocsolver_lwork = 0;
        CHECK_CUSOLVER(rocsolverDnZhegvd_bufferSize(
            rocsolver_handle,
            // TODO : handle
            ROCSOLVER_EIG_TYPE_1,
            ROCSOLVER_EIG_MODE_VECTOR,
            HIPBLAS_FILL_MODE_UPPER,
            nstart,
            hvec,
            ldh,
            sdum,
            ldh,
            e,
            &rocsolver_lwork));
        // cout<<"work_space: "<<rocsolver_lwork<<endl;
        hipblasDoubleComplex *rocsolver_work;
        CHECK_CUDA(hipMalloc((void**)&rocsolver_work, rocsolver_lwork*sizeof(hipblasDoubleComplex)));
        CHECK_CUSOLVER(rocsolverDnZhegvd(
            rocsolver_handle,
            ROCSOLVER_EIG_TYPE_1,
            ROCSOLVER_EIG_MODE_VECTOR,
            HIPBLAS_FILL_MODE_UPPER,
            nstart,
            hvec,
            ldh,
            sdum,
            ldh,
            e,
            rocsolver_work,
            rocsolver_lwork,
            device_info));
        CHECK_CUDA(hipFree(rocsolver_work));
        // CHECK_CUDA(hipFree(device_info));
    }
    else // nstart != nbands
    {
        // cout<<"not all"<<endl;
        CHECK_CUDA(hipMalloc((void**)&hdum, nstart*ldh*sizeof(hipblasDoubleComplex)));
        // hdum = hc;
        CHECK_CUDA(hipMemcpy(hdum, hc, nstart*ldh*sizeof(hipblasDoubleComplex), hipMemcpyDeviceToDevice));
        hipDeviceSynchronize();
        int rocsolver_lwork = 0;
        CHECK_CUSOLVER(rocsolver_zhegvdx_bufferSize(
            rocsolver_handle,
            ROCSOLVER_EIG_TYPE_1,
            ROCSOLVER_EIG_MODE_VECTOR,
            ROCSOLVER_EIG_RANGE_I,
            HIPBLAS_FILL_MODE_UPPER,
            nstart,
            hdum,
            ldh,
            sdum,
            ldh,
            0.0,
            0.0,
            1,
            nbands,
            &h_meig,
            e,
            &rocsolver_lwork));
        hipblasDoubleComplex *rocsolver_work;
        CHECK_CUDA(hipMalloc((void**)&rocsolver_work, rocsolver_lwork*sizeof(hipblasDoubleComplex)));
        CHECK_CUSOLVER(rocsolverDnZhegvdx(
            rocsolver_handle,
            rocblas_eform_ax,
            rocblas_evect_original,
            ROCSOLVER_EIG_RANGE_I,
            HIPBLAS_FILL_MODE_UPPER,
            nstart,
            hdum,
            ldh,
            sdum,
            ldh,
            0.0,
            0.0,
            1,
            nbands,
            &h_meig, // ?
            e,
            rocsolver_work,
            rocsolver_lwork,
            device_info
        ));
        CHECK_CUDA(hipFree(rocsolver_work));
    }
    hipDeviceSynchronize();
    // cout<<"end zhegv"<<endl;
    CHECK_CUDA(hipFree(device_info));
    // CHECK_CUDA(hipFree(rocsolver_work));
    // CHECK_CUSOLVER(rocsolverDnDestroy(rocsolver_handle));

	ModuleBase::timer::tick("Hamilt","diagH_CUSOLVER");
    return;
}

*/