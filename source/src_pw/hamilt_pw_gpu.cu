#include "tools.h"
#include "global.h"
#include "hamilt_pw_gpu.h"
#include "../module_base/blas_connector.h"
#include "../src_io/optical.h" // only get judgement to calculate optical matrix or not.
#include "myfunc.h"

__global__ void kernel_copy(int size, CUFFT_COMPLEX* dst, const CUFFT_COMPLEX *src)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x;
        dst[idx].y = src[idx].y;
    }
}

__global__ void kernel_get_tmhpsi(int size, CUFFT_COMPLEX *dst, const CUFFT_COMPLEX *src, double *g2kin)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x * g2kin[idx];
        dst[idx].y = src[idx].y * g2kin[idx];
    }
}

__global__ void kernel_add_tmhpsi(int size, CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src, int *index)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int p = index[idx];
    if(idx < size)
    {
        dst[idx].x += src[p].x;
        dst[idx].y += src[p].y;
    }
}

__global__ void kernel_addpp(CUFFT_COMPLEX *ps, double *deeq, const CUFFT_COMPLEX *becp, int nproj, int nprojx, int sum, int m, int nkb)
{
    int ip2 = blockDim.x * blockIdx.x + threadIdx.x;
    int ib = blockDim.y * blockIdx.y + threadIdx.y;
    if(ip2<nproj && ib<m)
    {
        ps[(sum+ip2) * m + ib].x = ps[(sum+ip2) * m + ib].y = 0;

        for(int ip=0; ip<nproj; ip++)
        {
            ps[(sum+ip2) * m + ib].x += deeq[ip * nprojx + ip2] * becp[ib * nkb + sum + ip].x;
            ps[(sum+ip2) * m + ib].y += deeq[ip * nprojx + ip2] * becp[ib * nkb + sum + ip].y;
        }
        // __syncthreads();
    }
}


int Hamilt_PW::moved = 0;

Hamilt_PW::Hamilt_PW()
{
    // hpsi = new complex<double>[1];
    // spsi = new complex<double>[1];
    // GR_index = new int[1];
    // Bec = new complex<double>[1];
    cudaMalloc((void**)&GR_index, sizeof(int)); // Only use this member now.
}

Hamilt_PW::~Hamilt_PW()
{
    // delete[] hpsi;
    // delete[] spsi;
    // delete[] GR_index;
    // delete[] Bec;
    cudaFree(GR_index);
}


void Hamilt_PW::allocate(
	const int &npwx,
	const int &npol,
	const int &nkb,
	const int &nrxx)
{
    TITLE("Hamilt_PW_GPU","allocate");

	assert(npwx > 0);
	assert(npol > 0);
	assert(nkb >=0);
	assert(nrxx > 0);

    // delete[] hpsi;
    // delete[] spsi;
    // delete[] GR_index;
    // delete[] Bec;

    cudaFree(GR_index);
    // this->hpsi = new complex<double> [npwx * npol];
    // this->spsi = new complex<double> [npwx * npol];
    cudaMalloc((void**)&GR_index, nrxx*sizeof(int));
    // this->Bec = new complex<double> [nkb];

    // ZEROS(this->hpsi, npwx * npol);
    // ZEROS(this->spsi, npwx * npol);
    // ZEROS(this->GR_index, nrxx);

    return;
}


void Hamilt_PW::init_k(const int ik)
{
    TITLE("Hamilt_PW_GPU","init_k");
	// mohan add 2010-09-30
	// (1) Which spin to use.
	if(GlobalV::NSPIN==2)
	{
		GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
	}

	// (2) Kinetic energy.
	GlobalC::wf.ekin(ik);

	// (3) Take the local potential.
	// cout<<"nrxx="<<GlobalC::pw.nrxx<<endl;

	for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);//mohan add 2007-11-12
	}

	// (4) Calculate nonlocal pseudopotential vkb
	//if (GlobalC::ppcell.nkb > 0 && !LINEAR_SCALING) xiaohui modify 2013-09-02
	if(GlobalC::ppcell.nkb > 0 && (GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")) //xiaohui add 2013-09-02. Attention...
	{
		GlobalC::ppcell.getvnl(ik);
	}

	// (5) The number of wave functions.
	GlobalC::wf.npw = GlobalC::kv.ngk[ik];

	// (6) The index of plane waves.
    int *GR_index_tmp = new int[GlobalC::pw.nrxx];
    for (int ig = 0;ig < GlobalC::wf.npw;ig++)
    {
        GR_index_tmp[ig] = GlobalC::pw.ig2fftw[ GlobalC::wf.igk(ik, ig) ];
    }
    // cout<<"init_K"<<endl;
    cudaMemcpy(this->GR_index, GR_index_tmp, GlobalC::pw.nrxx*sizeof(int), cudaMemcpyHostToDevice);
    delete [] GR_index_tmp;

    // (7) ik
	GlobalV::CURRENT_K = ik;

    return;
}


//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
void Hamilt_PW::diagH_subspace(
    const int ik,
    const int nstart,
    const int n_band,
    const ComplexMatrix &psi,
    ComplexMatrix &evc,
    double *en)
{
    TITLE("Hamilt_PW","diagH_subspace");
    timer::tick("Hamilt_PW","diagH_subspace");

	assert(nstart!=0);
	assert(n_band!=0);

    ComplexMatrix hc(nstart, nstart);
    ComplexMatrix sc(nstart, nstart);
    ComplexMatrix hvec(nstart,n_band);

	int dmin=0;
	int dmax=0;
	const int npw = GlobalC::kv.ngk[ik];

	if(GlobalV::NSPIN != 4)
	{
		dmin= npw;
		dmax = GlobalC::wf.npwx;
	}
	else
	{
		dmin = GlobalC::wf.npwx*GlobalV::NPOL;
		dmax = GlobalC::wf.npwx*GlobalV::NPOL;
	}

	//qianrui improve this part 2021-3-14
	std::complex<double> *aux=new std::complex<double> [dmax*nstart];
	std::complex<double> *paux = aux;
	std::complex<double> *ppsi = psi.c;

	//qianrui replace it
	this->h_psi(psi.c, aux, nstart);

	char trans1 = 'C';
	char trans2 = 'N';
	zgemm_(&trans1,&trans2,&nstart,&nstart,&dmin,&ONE,psi.c,&dmax,aux,&dmax,&ZERO,hc.c,&nstart);
	hc=transpose(hc,false);

	zgemm_(&trans1,&trans2,&nstart,&nstart,&dmin,&ONE,psi.c,&dmax,psi.c,&dmax,&ZERO,sc.c,&nstart);
	sc=transpose(sc,false);

	delete []aux;

	// Peize Lin add 2019-03-09
#ifdef __LCAO
	if(GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		auto add_Hexx = [&](const double alpha)
		{
			for (int m=0; m<nstart; ++m)
			{
				for (int n=0; n<nstart; ++n)
				{
					hc(m,n) += alpha * GlobalC::exx_lip.get_exx_matrix()[ik][m][n];
				}
			}
		};
		if( 5==GlobalC::xcf.iexch_now && 0==GlobalC::xcf.igcx_now )				// HF
		{
			add_Hexx(1);
		}
		else if( 6==GlobalC::xcf.iexch_now && 8==GlobalC::xcf.igcx_now )			// PBE0
		{
			add_Hexx(GlobalC::exx_global.info.hybrid_alpha);
		}
		else if( 9==GlobalC::xcf.iexch_now && 12==GlobalC::xcf.igcx_now )			// HSE
		{
			add_Hexx(GlobalC::exx_global.info.hybrid_alpha);
		}
	}
#endif

	if(GlobalV::NPROC_IN_POOL>1)
	{
		Parallel_Reduce::reduce_complex_double_pool( hc.c, nstart*nstart );
		Parallel_Reduce::reduce_complex_double_pool( sc.c, nstart*nstart );
	}

	// after generation of H and S matrix, diag them
    GlobalC::hm.diagH_LAPACK(nstart, n_band, hc, sc, nstart, en, hvec);


	// Peize Lin add 2019-03-09
#ifdef __LCAO
	if("lcao_in_pw"==GlobalV::BASIS_TYPE)
	{
		switch(GlobalC::exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_lip.k_pack->hvec_array[ik] = hvec;
				break;
		}
	}
#endif

    //=======================
    //diagonize the H-matrix
    //=======================

// for tests
/*
		std::cout << std::setprecision(3);
		out.printV3(GlobalV::ofs_running,GlobalC::kv.kvec_c[ik]);
		out.printcm_norm("sc",sc,1.0e-4);
		out.printcm_norm("hvec",hvec,1.0e-4);
		out.printcm_norm("hc",hc,1.0e-4);
		std::cout << std::endl;
*/

	std::cout << std::setprecision(5);

//--------------------------
// KEEP THIS BLOCK FOR TESTS
//--------------------------
/*
	std::cout << "  hc matrix" << std::endl;
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			double a = hc(i,j).real();
			if(abs(a) < 1.0e-5) a = 0;
			std::cout << std::setw(6) << a;
		}
		std::cout << std::endl;
	}

	std::cout << "  sc matrix" << std::endl;
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			double a = sc(i,j).real();
			if(abs(a) < 1.0e-5) a = 0;
			std::cout << std::setw(6) << a;
		}
		std::cout << std::endl;
	}

	std::cout << "\n Band Energy" << std::endl;
	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		std::cout << " e[" << i+1 << "]=" << en[i] * Ry_to_eV << std::endl;
	}
*/
//--------------------------
// KEEP THIS BLOCK FOR TESTS
//--------------------------


	if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") && GlobalV::CALCULATION=="nscf" && !Optical::opt_epsilon2)
	{
		GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
	}
	else if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		&& ( GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax")) //pengfei 2014-10-13
	{
		// because psi and evc are different here,
		// I think if psi and evc are the same,
		// there may be problems, mohan 2011-01-01
		char transa = 'N';
		char transb = 'T';
		zgemm_( &transa,
				&transb,
				&dmax, // m: row of A,C
				&n_band, // n: col of B,C
				&nstart, // k: col of A, row of B
				&ONE, // alpha
				psi.c, // A
				&dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
				hvec.c, // B
				&n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
				&ZERO,  // belta
				evc.c, // C
				&dmax ); // LDC: if(N) max(1, m)
	}
	else
	{
		// As the evc and psi may refer to the same matrix, we first
		// create a temporary matrix to story the result. (by wangjp)
		// qianrui improve this part 2021-3-13
		char transa = 'N';
		char transb = 'T';
		ComplexMatrix evctmp(n_band, dmin,false);
		zgemm_(&transa,&transb,&dmin,&n_band,&nstart,&ONE,psi.c,&dmax,hvec.c,&n_band,&ZERO,evctmp.c,&dmin);
		for(int ib=0; ib<n_band; ib++)
		{
			for(int ig=0; ig<dmin; ig++)
			{
				evc(ib,ig) = evctmp(ib,ig);
			}
		}
	}
    //out.printr1_d("en",en,n_band);

//	std::cout << "\n bands" << std::endl;
//	for(int ib=0; ib<n_band; ib++)
//	{
//		std::cout << " ib=" << ib << " " << en[ib] * Ry_to_eV << std::endl;
//	}

    //out.printcm_norm("hvec",hvec,1.0e-8);

    timer::tick("Hamilt_PW","diagH_subspace");
    return;
}

void Hamilt_PW::h_1psi( const int npw_in, const CUFFT_COMPLEX *psi,
                        CUFFT_COMPLEX *hpsi, CUFFT_COMPLEX *spsi)
{
    this->h_psi(psi, hpsi);

    int thread = 512;
    int block = npw_in / thread + 1;
    kernel_copy<<<thread, block>>>(npw_in, spsi, psi);
    return;
}

void Hamilt_PW::s_1psi(const int dim, const CUFFT_COMPLEX *psi, CUFFT_COMPLEX *spsi)
{
    cudaMemcpy(spsi, psi, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
    return;
}

void Hamilt_PW::h_psi(const CUFFT_COMPLEX *psi_in, CUFFT_COMPLEX *hpsi, const int m)
{
    timer::tick("Hamilt_PW_GPU","h_psi");
    // int i = 0;
    // int j = 0;
    // int ig= 0;

	//if(NSPIN!=4) ZEROS(hpsi, wf.npw);
	//else ZEROS(hpsi, wf.npwx * NPOL);//added by zhengdy-soc
	int dmax = GlobalC::wf.npwx * GlobalV::NPOL;

    // cout<<"dim inside="<<GlobalC::wf.npwx * GlobalV::NPOL<<endl;

	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	CUFFT_COMPLEX *tmhpsi;
	const CUFFT_COMPLEX *tmpsi_in;
    timer::tick("Hamilt_PW_GPU","kinetic");
 	if(GlobalV::T_IN_H)
	{
        tmhpsi = hpsi;
        tmpsi_in = psi_in;

        double* d_g2kin;
        cudaMalloc((void**)&d_g2kin, GlobalC::wf.npwx*sizeof(double));
        cudaMemcpy(d_g2kin, GlobalC::wf.g2kin, GlobalC::wf.npw*sizeof(double), cudaMemcpyHostToDevice);
        for(int ib = 0 ; ib < m; ++ib)
        {
            // cout<<"in hpsi-Kinetic, iband = "<<ib<<endl;

            int thread = 512;
            int block = GlobalC::wf.npw / thread + 1;
            kernel_get_tmhpsi<<<block, thread>>>(GlobalC::wf.npw, tmhpsi, tmpsi_in, d_g2kin);

            // if(GlobalC::::NSPIN==4){
            //     for(ig=GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
            //     {
            //         tmhpsi[ig] = 0;
            //     }
            //     tmhpsi += GlobalC::wf.npwx;
            //     tmpsi_in += GlobalC::wf.npwx;
            //     for (ig = 0;ig < GlobalC::wf.npw ;++ig)
            //     {
            //         tmhpsi[ig] = GlobalC::wf.g2kin[ig] * tmpsi_in[ig];
            //     }
            //     // TODO: setup with 0
            //     for(ig=GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
            //     {
            //         tmhpsi[ig] =0;
            //     }
            // }

            tmhpsi += GlobalC::wf.npwx;
            tmpsi_in += GlobalC::wf.npwx;
        }
        cudaFree(d_g2kin);
	}
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
    }
    timer::tick("Hamilt_PW_GPU","kinetic");

	//------------------------------------
	//(2) the local potential.
	//-----------------------------------
	timer::tick("Hamilt_PW_GPU","vloc");
    //  ...
	if(GlobalV::VL_IN_H)
	{
        tmhpsi = hpsi;
        tmpsi_in = psi_in;
        // int *d_GR_index;
        double *d_vr_eff1;
        CUFFT_COMPLEX *d_porter;

        // cudaMalloc((void**)&d_GR_index, GlobalC::wf.npwx * sizeof(int));
        cudaMalloc((void**)&d_vr_eff1, GlobalC::pw.nrxx * sizeof(double));
        cudaMalloc((void**)&d_porter, GlobalC::pw.nrxx * sizeof(CUFFT_COMPLEX));

        cudaMemcpy(d_vr_eff1, GlobalC::pot.vr_eff1, GlobalC::pw.nrxx*sizeof(double), cudaMemcpyHostToDevice);
        // cout<<"NSPIN = "<<GlobalV::NSPIN<<endl;
        for(int ib = 0 ; ib < m; ++ib)
        {
            // cout<<"in hpsi:loacl_pot, iband = "<<ib<<endl;
            // if(NSPIN!=4){
            // ZEROS( UFFT.porter, pw.nrxx);
            cudaMemset(d_porter, 0, GlobalC::pw.nrxx * sizeof(CUFFT_COMPLEX));

            GlobalC::UFFT.RoundTrip_GPU( tmpsi_in, d_vr_eff1, GR_index, d_porter );

            // for (j = 0;j < wf.npw;j++)
            // {
            //     tmhpsi[j] += UFFT.porter[ GR_index[j] ];
            // }
            int thread = 512;
            int block = GlobalC::wf.npw / thread + 1;
            kernel_add_tmhpsi<<<block, thread>>>(GlobalC::wf.npw, tmhpsi, d_porter, GR_index);

            tmhpsi += dmax;
            tmpsi_in += dmax;
        }
        // cudaFree(d_GR_index);
        cudaFree(d_vr_eff1);
        cudaFree(d_porter);
	}
	timer::tick("Hamilt_PW_GPU","vloc");
	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	timer::tick("Hamilt_PW_GPU","vnl");
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
    }

    if(GlobalV::VNL_IN_H)
	{
        if ( GlobalC::ppcell.nkb > 0)
        {
            int nkb = GlobalC::ppcell.nkb;
            CUFFT_COMPLEX *becp;
            CUFFT_COMPLEX *d_vkb_c;
            cudaMalloc((void**)&becp, GlobalV::NPOL*m*nkb*sizeof(CUFFT_COMPLEX));
            cudaMalloc((void**)&d_vkb_c, GlobalC::wf.npwx*nkb*sizeof(CUFFT_COMPLEX));

            cudaMemcpy(d_vkb_c, GlobalC::ppcell.vkb.c, GlobalC::wf.npwx*nkb*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);

            cublasOperation_t transa = CUBLAS_OP_C;
            cublasOperation_t transb = CUBLAS_OP_N;
            cublasHandle_t handle;
            cublasCreate(&handle);

            CUFFT_COMPLEX ONE, ZERO;
            ONE.y = ZERO.x = ZERO.y = 0.0;
            ONE.x = 1.0;
            // NEG_ONE.x = -1.0;

            if(m==1 && GlobalV::NPOL==1)
            {
                int inc = 1;
                cublasZgemv(handle, transa, GlobalC::wf.npw, nkb, &ONE, d_vkb_c, GlobalC::wf.npwx, psi_in, inc, &ZERO, becp, inc);

            }
            else
            {
                int npm = GlobalV::NPOL * m;
                cublasZgemm(handle, transa, transb, nkb, npm, GlobalC::wf.npw, &ONE, d_vkb_c, GlobalC::wf.npwx, psi_in, GlobalC::wf.npwx, &ZERO, becp, nkb);

            }

            // complex<double> *hpsi_cpu = new complex<double>[GlobalC::wf.npw*GlobalV::NPOL];
            // complex<double> *becp_cpu = new complex<double>[GlobalV::NPOL*m*nkb];

            // cudaMemcpy(becp_cpu, becp, GlobalV::NPOL*m*nkb*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);

            // cudaMemcpy(hpsi_cpu, hpsi, GlobalC::wf.npw*GlobalV::NPOL*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);

            // this->add_nonlocal_pp(hpsi_cpu, becp_cpu, m);

            // cudaMemcpy(hpsi, hpsi_cpu, GlobalC::wf.npw*GlobalV::NPOL*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);

            // delete [] hpsi_cpu;
            // delete [] becp_cpu;

            this->add_nonlocal_pp(hpsi, becp, d_vkb_c, m);

            cublasDestroy(handle);
            cudaFree(becp);
            cudaFree(d_vkb_c);
            // cout<<"nonlocal end"<<endl;

        }
    }

    timer::tick("Hamilt_PW_GPU","vnl");

    //------------------------------------
	// (4) the metaGGA part
	//------------------------------------
    // TODO: add metaGGA part

    timer::tick("Hamilt_PW_GPU","h_psi");
    return;
}


void Hamilt_PW::add_nonlocal_pp(
	CUFFT_COMPLEX *hpsi_in,
	const CUFFT_COMPLEX *becp,
    const CUFFT_COMPLEX *d_vkb_c,
	const int m)
{
    timer::tick("Hamilt_PW_GPU","add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	// complex<double> *ps  = new complex<double> [nkb * GlobalV::NPOL * m];
    // ZEROS(ps, GlobalV::NPOL * m * nkb);
    CUFFT_COMPLEX *ps;
    cudaMalloc((void**)&ps, nkb * GlobalV::NPOL * m * sizeof(CUFFT_COMPLEX));
    cudaMemset(ps, 0, GlobalV::NPOL * m * sizeof(CUFFT_COMPLEX));

    int sum = 0;
    int iat = 0;
    // if(GlobalV::NSPIN!=4)
	// {
    for (int it=0; it<GlobalC::ucell.ntype; it++)
    {

        const int nproj = GlobalC::ucell.atoms[it].nh;
        const int nprojx = GlobalC::ppcell.nhm;
        double *cur_deeq;
        cudaMalloc((void**)&cur_deeq, nprojx*nprojx*sizeof(double));
        for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
        {
            cudaMemcpy(cur_deeq, &(GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, 0, 0)),
                nprojx*nprojx*sizeof(double), cudaMemcpyHostToDevice);

            int thread_x = 16;
            dim3 thread(thread_x, thread_x);
            dim3 block((nproj+thread_x-1)/thread_x, (m+thread_x-1)/thread_x);
            // dim3 block(1, 1, 1);

            kernel_addpp<<<block, thread>>>(ps, cur_deeq, becp, nproj, nprojx, sum, m, nkb);

            sum += nproj;
            ++iat;
        } //end na
    } //end nt
	// }

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
    }

    cublasOperation_t transa = CUBLAS_OP_N;
    cublasOperation_t transb = CUBLAS_OP_T;
    cublasHandle_t handle;
    cublasCreate(&handle);
    CUFFT_COMPLEX ONE;
    ONE.y = 0.0;
    ONE.x = 1.0;
	if(GlobalV::NPOL==1 && m==1)
	{
		int inc = 1;
        cublasZgemv(handle,
            transa,
            GlobalC::wf.npw,
            GlobalC::ppcell.nkb,
            &ONE,
            d_vkb_c,
            GlobalC::wf.npwx,
			ps,
			inc,
			&ONE,
			hpsi_in,
			inc);
	}
	else
	{
		int npm = GlobalV::NPOL*m;
        cublasZgemm(handle,
            transa,
            transb,
            GlobalC::wf.npw,
            npm,
            GlobalC::ppcell.nkb,
            &ONE,
            d_vkb_c,
            GlobalC::wf.npwx,
            ps,
            npm,
            &ONE,
            hpsi_in,
            GlobalC::wf.npwx);
	}

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
    }

	// delete[] ps;
    cudaFree(ps);
    timer::tick("Hamilt_PW_GPU","add_nonlocal_pp");
    return;
}
