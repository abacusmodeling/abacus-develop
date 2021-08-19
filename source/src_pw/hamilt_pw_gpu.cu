#include "tools.h"
#include "global.h"
#include "hamilt_pw.h"
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


int Hamilt_PW_GPU::moved = 0;

Hamilt_PW_GPU::Hamilt_PW_GPU()
{
    // hpsi = new complex<double>[1];
    // spsi = new complex<double>[1];
    // GR_index = new int[1];
    // Bec = new complex<double>[1];
    cudaMalloc((void**)&GR_index, sizeof(int)); // Only use this member now.
}

Hamilt_PW_GPU::~Hamilt_PW_GPU()
{
    // delete[] hpsi;
    // delete[] spsi;
    // delete[] GR_index;
    // delete[] Bec;
    cudaFree(GR_index);
}


void Hamilt_PW_GPU::allocate(
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


void Hamilt_PW_GPU::init_k(const int ik)
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
    return;
}

void Hamilt_PW_GPU::s_1psi(const int dim, const CUFFT_COMPLEX *psi, CUFFT_COMPLEX *spsi)
{
    cudaMemcpy(spsi, psi, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
    return;
}

void Hamilt_PW_GPU::h_1psi( const int npw_in, const CUFFT_COMPLEX *psi,
            CUFFT_COMPLEX *hpsi, CUFFT_COMPLEX *spsi)
{
    this->h_psi(psi, hpsi);

    int thread = 512;
    int block = npw_in / thread + 1;
    kernel_copy<<<thread, block>>>(npw_in, spsi, psi);
    return;
}

void Hamilt_PW_GPU::h_psi(const CUFFT_COMPLEX *psi_in, CUFFT_COMPLEX *hpsi, const int m)
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

            this->add_nonlocal_pp_gpu(hpsi, becp, d_vkb_c, m);

            cublasDestroy(handle);
            cudaFree(becp);
            cudaFree(d_vkb_c);
            // cout<<"nonlocal end"<<endl;

        }
    }
    
    timer::tick("Hamilt_PW_GPU","vnl");
    timer::tick("Hamilt_PW_GPU","h_psi");
    return;
}

void Hamilt_PW_GPU::add_nonlocal_pp(
	complex<double> *hpsi_in,
	const complex<double> *becp,
	const int m)
{
    timer::tick("Hamilt_PW_GPU","add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	complex<double> *ps  = new complex<double> [nkb * GlobalV::NPOL * m];
    ZEROS(ps, GlobalV::NPOL * m * nkb);

    int sum = 0;
    int iat = 0;
    if(GlobalV::NSPIN!=4)
	{
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<nproj; ip++)
				{
                    for (int ip2=0; ip2<nproj; ip2++)
                    {
                        for(int ib = 0; ib < m ; ++ib)
                        {
                            ps[(sum + ip2) * m + ib] += 
                            GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2) 
                            * becp[ib * nkb + sum + ip];
                        }//end ib
                    }// end ih
				}//end jh 
				sum += nproj;
				++iat;
			} //end na
		} //end nt
	}
	else
	{
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			int psind=0;
			int becpind=0;
			complex<double> becp1=complex<double>(0.0,0.0);
			complex<double> becp2=complex<double>(0.0,0.0);

			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.

				for (int ip=0; ip<nproj; ip++)
				{
					for (int ip2=0; ip2<nproj; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							psind = (sum+ip2) * 2 * m + ib * 2;
							becpind = ib*nkb*2 + sum + ip;
							becp1 =  becp[becpind];
							becp2 =  becp[becpind + nkb];
							ps[psind] += GlobalC::ppcell.deeq_nc(0, iat, ip2, ip) * becp1
								+GlobalC::ppcell.deeq_nc(1, iat, ip2, ip) * becp2;
							ps[psind +1] += GlobalC::ppcell.deeq_nc(2, iat, ip2, ip) * becp1
								+GlobalC::ppcell.deeq_nc(3, iat, ip2, ip) * becp2;
						}//end ib
					}// end ih
				}//end jh
				sum += nproj;
				++iat;
			} //end na
		} //end nt
	}

	/*
    for (int ig=0;ig<GlobalC::wf.npw;ig++)
    {
        for (int i=0;i< GlobalC::ppcell.nkb;i++)
        {
            hpsi_in[ig]+=ps[i]*GlobalC::ppcell.vkb(i,ig);
        }
    }
	*/

	// use simple method.
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//qianrui optimize 2021-3-31
	char transa = 'N';
	char transb = 'T';
	if(GlobalV::NPOL==1 && m==1)
	{
		int inc = 1;
		zgemv_(&transa, 
			&GlobalC::wf.npw, 
			&GlobalC::ppcell.nkb, 
			&ONE, 
			GlobalC::ppcell.vkb.c, 
			&GlobalC::wf.npwx, 
			ps, 
			&inc, 
			&ONE, 
			hpsi_in, 
			&inc);
	}
	else
	{
		int npm = GlobalV::NPOL*m;
		zgemm_(&transa,
			&transb,
			&GlobalC::wf.npw,
			&npm,
			&GlobalC::ppcell.nkb,
			&ONE,
			GlobalC::ppcell.vkb.c,
			&GlobalC::wf.npwx,
			ps,
			&npm,
			&ONE,
			hpsi_in,
			&GlobalC::wf.npwx);
	}

	delete[] ps;
    timer::tick("Hamilt_PW_GPU","add_nonlocal_pp");
    return;
}

void Hamilt_PW_GPU::add_nonlocal_pp_gpu(
	CUFFT_COMPLEX *hpsi_in,
	const CUFFT_COMPLEX *becp,
    const CUFFT_COMPLEX *d_vkb_c,
	const int m)
{
    timer::tick("Hamilt_PW_GPU","add_nonlocal_pp_gpu");

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
    timer::tick("Hamilt_PW_GPU","add_nonlocal_pp_gpu");
    return;
}
