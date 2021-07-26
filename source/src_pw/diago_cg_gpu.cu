#include "diago_cg_gpu.h"
#include "cuda_runtime.h"
#include "../src_pw/global.h"

int Diago_CG_GPU::moved = 0;

__global__ void kernel_normalization(CUFFT_COMPLEX *data, int size, double norm)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        data[idx].x /= norm;
        data[idx].y /= norm;
    }
}

__global__ void kernel_precondition(CUFFT_COMPLEX *res, const CUFFT_COMPLEX *data, const int size, const double *P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        res[idx].x = data[idx].x / P[idx];
        res[idx].y = data[idx].y / P[idx];
    }
}

__global__ void kernel_precondition_inverse(CUFFT_COMPLEX *res, const CUFFT_COMPLEX *data, const int size, const double *P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        res[idx].x = data[idx].x * P[idx];
        res[idx].y = data[idx].y * P[idx];
    }
}

__global__ void kernel_get_gredient(CUFFT_COMPLEX *g, CUFFT_COMPLEX *ppsi, int size, double lambda)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        g[idx].x -= lambda * ppsi[idx].x;
        g[idx].y -= lambda * ppsi[idx].y;
    }
}

__global__ void kernel_get_gammacg(int size, CUFFT_COMPLEX *dst, const CUFFT_COMPLEX *src, double gamma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = dst[idx].x * gamma + src[idx].x;
        dst[idx].y = dst[idx].y * gamma + src[idx].y;
    }
}

__global__ void kernel_get_normacg(int size, CUFFT_COMPLEX *dst, const CUFFT_COMPLEX *src, double norma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = dst[idx].x - norma * src[idx].x;
        dst[idx].y = dst[idx].y - norma * src[idx].y;
    }
}

__global__ void kernel_multi_add(CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src1, double a1, const CUFFT_COMPLEX *src2, double a2, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src1[idx].x * a1 + src2[idx].x * a2;
        dst[idx].y = src1[idx].y * a1 + src2[idx].y * a2;
    }
}

// __global__ void fake_s_1psi(int size, const CUFFT_COMPLEX *src, CUFFT_COMPLEX *dst)
// {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if(idx < size)
//     {
//         dst[idx].x = src[idx].x;
//         dst[idx].y = src[idx].y;
//     }
// }

Diago_CG_GPU::Diago_CG_GPU()
{
    test_cg=0;
}

Diago_CG_GPU::~Diago_CG_GPU() {}

void Diago_CG_GPU::diag
(
    CUFFT_COMPLEX *phi, // matrix nband*dim 
    double *e, //
    const int &dim,//
    const int &dmx,//
    const int &n_band,//能带数目
    const double *precondition,//
    const double &eps,//
    const int &maxter,//最大cg法迭代次数
    const bool &reorder,//是否排序
    int &notconv,//
    double &avg_iter//
)
{
    if (test_cg==1) TITLE("Diago_CG_GPU","ccgdiagg");
    timer::tick("Diago_CG_GPU","diag");

    avg_iter = 0.0;
    notconv = 0;
    // ZEROS(e, n_band);

    //-------------------------------------------------------------------
    // "poor man" iterative diagonalization of a complex hermitian matrix
    // through preconditioned conjugate gradient algorithm
    // Band-by-band algorithm with minimal use of memory
    // Calls h_1phi and s_1phi to calculate H|phi> and S|phi>
    // Works for generalized eigenvalue problem (US pseudopotentials) as well
    //-------------------------------------------------------------------
    
    CUFFT_COMPLEX *sphi;
    CUFFT_COMPLEX *scg;
    CUFFT_COMPLEX *hphi;
    CUFFT_COMPLEX *g;
    CUFFT_COMPLEX *cg;
    CUFFT_COMPLEX *g0;
    CUFFT_COMPLEX *pphi;
    CUFFT_COMPLEX *lagrange;
    CUFFT_COMPLEX *phi_m;

    cout << "Hello, CG!" << endl;
    cout << "CG Dim = " << dim << " & " << dmx << endl;

    CUFFT_COMPLEX *test_malloc;
    cudaMalloc((void**)&test_malloc, dim*sizeof(CUFFT_COMPLEX));

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
    cudaFree(test_malloc);

    cudaMalloc((void**)&sphi, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&scg, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&hphi, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&g, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&cg, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&g0, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&pphi, dim * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&lagrange, n_band * sizeof(CUFFT_COMPLEX));
    cudaMalloc((void**)&phi_m, dim * sizeof(CUFFT_COMPLEX));

    err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

	// Init with ZERO ...
    // double em_host = 0;

    for (int m=0; m<n_band; m++) 
    {
        // if (test_cg>2) ofs_running << "Diagonal Band : " << m << endl;
        cout << "Diagonal Band : " << m <<" of "<<n_band<< endl;

        // cout<<"!!!!!!!!!!!!  Test begin FOR iter ... "<<endl;
        // complex<double> *test_for = new complex<double>[15];
        // cudaMemcpy(test_for, &phi[(m+1)*dim], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
        // for(int i=0;i<15;i++)
        // {
        //     cout<<test_for[i].real()<<" "<<test_for[i].imag()<<endl;
        // }
        // delete [] test_for;
        // for (int i=0; i<dim; i++) phi_m[i] = phi[m*dim + i];
        
        // cout<<"Test Phi_m ... "<<endl;
        // complex<double> *test_phim = new complex<double>[15];
        // cudaMemcpy(test_phim, &phi[m*dim], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
        // for(int i=0;i<15;i++)
        // {
        //     cout<<test_phim[i].real()<<" "<<test_phim[i].imag()<<endl;
        // }
        // delete [] test_phim;

        cudaMemcpy(phi_m, &phi[m*dmx], dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);

        err = cudaGetLastError();
        if (err != cudaSuccess)
            cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

        // pw.s_1psi(dim, phi_m, sphi); // TODO
        // pw难整合，先用copy替代，恒等变换 
        cudaMemcpy(sphi, phi_m, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);

        err = cudaGetLastError();
        if (err != cudaSuccess)
            cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
        
        this->schmit_orth(dim, dmx, m, phi, sphi, phi_m);
        // pw.h_1psi(dim, phi_m, hphi, sphi);

        // To CPU
        // 涉及到complex模板类和cufftComplex的转换
        complex<double> *phi_m_cpu = new complex<double>[dim];
        complex<double> *hphi_cpu = new complex<double>[dim];
        complex<double> *sphi_cpu = new complex<double>[dim];
        
        cudaMemcpy(phi_m_cpu, phi_m, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost); // TODO
        cudaMemcpy(hphi_cpu, hphi, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
        cudaMemcpy(sphi_cpu, sphi, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);

        // cout<<"Info before h_1psi result: ...."<<endl;
        // for(int i=0;i<10;i++)
        // {
        //     cout<<phi_m_cpu[i].real()<<" "<<phi_m_cpu[i].imag()<<endl;
        // }

        GlobalC::hm.hpw.h_1psi(dim, phi_m_cpu, hphi_cpu, sphi_cpu);

        // cout<<"Info of h_1psi result: ...."<<endl;
        // for(int i=0;i<10;i++)
        // {
        //     cout<<hphi_cpu[i].real()<<" "<<hphi_cpu[i].imag()<<endl;
        // }
        
        // To GPU
        cudaMemcpy(phi_m, phi_m_cpu, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);
        cudaMemcpy(hphi, hphi_cpu, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);
        cudaMemcpy(sphi, sphi_cpu, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);

        err = cudaGetLastError();
        if (err != cudaSuccess)
            cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

        delete [] phi_m_cpu;
        delete [] hphi_cpu;
        delete [] sphi_cpu;

        cout<<"hpsi end"<<endl;

        cout<<"before ddot"<<endl;

        double em_host = 0;
        cout<<"EM_host1: "<<em_host<<endl;
        em_host = ddot_real(dim, phi_m, hphi);
        cout<<"EM_host2: "<<em_host<<endl;

        // cout<<"ddot!"<<endl;
        // cudaMemcpyToSymbol(&e[m], &em_host, sizeof(double));
        cudaMemcpy(&e[m], &em_host, sizeof(double), cudaMemcpyHostToDevice);
        cout<<"ddot real end"<<endl;

        int iter = 0;
        double gg_last = 0.0;
        double cg_norm = 0.0;
        double theta = 0.0;
        bool converged = false;
        // cg iteration

        for (iter = 0;iter < maxter; iter++)
        {
            cout<<"********iter*******:"<<iter<<endl;
            this->calculate_gradient( precondition, dim, hphi, sphi, g, pphi );
            cout<<"cal grad end"<<endl;

            this->orthogonal_gradient( dim, dmx, g, scg, lagrange, phi, m );
            cout<<"orth grad end"<<endl;
            this->calculate_gamma_cg( iter, dim, precondition, g, scg, 
			    g0, cg, gg_last, cg_norm, theta, phi_m);// scg used as sg
            cout<<"calcu gama end"<<endl;
            converged = this->update_psi( dim, cg_norm, theta, pphi, cg, scg, phi_m , 
			    em_host, eps, hphi, sphi); // pphi is used as hcg
            cudaMemcpy(&e[m], &em_host, sizeof(double), cudaMemcpyHostToDevice);
            cout<<"update psi end"<<endl;
            if ( converged ) break;
        }//end iter

        // for (int i = 0;i < dim;i++)
        // {
        //     phi[m * dim + i] = phi_m[i];
        // }
        
        cudaMemcpy(&phi[m*dmx], phi_m, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
        
        //经过maxter之后仍未收敛 notconv++

        // cout<<"****************  Test m+1 iter ... ***********"<<endl;
        // complex<double> *test_iter = new complex<double>[15];
        // cudaMemcpy(test_iter, &phi[(m+1)*dim], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
        // for(int i=0;i<15;i++)
        // {
        //     cout<<test_iter[i].real()<<" "<<test_iter[i].imag()<<endl;
        // }
        // delete [] test_iter;

        if (!converged)
        {
            ++notconv;
        }

        avg_iter += static_cast<double>(iter) + 1.00;
        // cout<<"iter:"<<iter<<endl;

        // reorder eigenvalues if they are not in the right order
        // (this CAN and WILL happen in not-so-special cases)

        if (m > 0 && reorder)
        {
			NOTE("reorder bands!");
            cout<<"reorder bands ..."<<endl;
            // TODO
            double* e_host;
            e_host = (double*)malloc(n_band*sizeof(double));
            ZEROS(e_host, n_band);
            cudaMemcpy(e_host, e, n_band*sizeof(double), cudaMemcpyDeviceToHost);
            
            // cout<<"e info before reorder .. "<<endl;
            // for(int i=0;i<10;i++)
            // {
            //     cout<<e_host[i]<<endl;
            // }
            
            if (e_host[m]-e_host[m-1]<-2.0*eps)
            {
                // if the last calculated eigenvalue is not the largest...
                int i=0;
                for (i=m-2; i>= 0; i--)
                {
                    if (e_host[m]-e_host[i]>2.0*eps) break;
                }
                i++;
                moved++;

                // last calculated eigenvalue should be in the i-th position: reorder
                double e0 = e_host[m];
                //dcopy(phi, m, pphi);

                // for (int ig=0;ig<dim;ig++)
                // {
                //     // pphi[ig]=phi(m,ig);
                //     pphi[ig] = phi[m*dim+ig];
                // }
                cudaMemcpy(pphi, &phi[m*dmx], dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);

                for (int j = m;j >= i + 1;j--)
                {
                    e_host[j]=e_host[j-1];
                    // for (int ig=0;ig<dim;ig++)
                    // {
                    //     // phi(j,ig) = phi(j-1,ig);
                    //     phi[j*dim+ig] = phi[(j-1)*dim+ig];
                    // }
                    cudaMemcpy(&phi[j*dmx], &phi[(j-1)*dmx], dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
                }

                e_host[i] = e0;
                //dcopy(pphi, phi, i);
                // for (int ig=0;ig<dim;ig++)
                // {
                //     // phi(i,ig) = pphi[ig];
                //     phi[i*dim+ig] = pphi[ig];
                // }
                cudaMemcpy(&phi[i*dmx], pphi, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
                // this procedure should be good if only a few inversions occur,
                // extremely inefficient if eigenvectors are often in bad order
                // (but this should not happen)
            } // endif

            cudaMemcpy(e, e_host, n_band*sizeof(double), cudaMemcpyHostToDevice);
            delete [] e_host;
        } //end reorder

    }//end m

    avg_iter /= n_band;

    // delete [] lagrange;
    // delete [] pphi;
    // delete [] g0;
    // delete [] cg;
    // delete [] g;
    // delete [] hphi;
    // delete [] scg;
    // delete [] sphi;
    // delete [] phi_m;

    cudaFree(lagrange);
    cudaFree(pphi);
    cudaFree(g0);
    cudaFree(cg);
    cudaFree(g);
    cudaFree(hphi);
    cudaFree(scg);
    cudaFree(sphi);
    cudaFree(phi_m);

    timer::tick("Diago_CG_GPU","diag");
    return;
} // end subroutine ccgdiagg


void Diago_CG_GPU::calculate_gradient(
    const double* precondition, const int dim,
    const CUFFT_COMPLEX *hpsi, const CUFFT_COMPLEX *spsi,
    CUFFT_COMPLEX *g, CUFFT_COMPLEX *ppsi)
{
    if (test_cg==1) TITLE("Diago_CG_GPU","calculate_gradient");
    //timer::tick("Diago_CG_GPU","grad");

    int thread = 512;
    int block = dim / thread + 1;

    // kernel_precondition(data, res, size, precondition)
    // (2) PH|psi> : g[i] = hpsi[i]/precondition[i]
    kernel_precondition<<<block, thread>>>(g, hpsi, dim, precondition);
    // (3) PS|psi> : ppsi[i] = spsi[i]/precondition[i]
    kernel_precondition<<<block, thread>>>(ppsi, spsi, dim, precondition);

    // Update lambda !
    // (4) <psi|SPH|psi >
    const double eh = this->ddot_real(dim, spsi, g);
    // (5) <psi|SPS|psi >
    const double es = this->ddot_real(dim, spsi, ppsi);
    const double lambda = eh / es;
    
    // Update g !
    kernel_get_gredient<<<block, thread>>>(g, ppsi, dim, lambda);
    // kernel_multi_add<<<block, thread>>>(g, g, 1, ppsi, -lambda, dim);
    //timer::tick("Diago_CG_GPU","grad");
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;
    return;
}


void Diago_CG_GPU::orthogonal_gradient( const int &dim, const int &dmx,
                                    CUFFT_COMPLEX *g, CUFFT_COMPLEX *sg, CUFFT_COMPLEX *lagrange,
                                    const CUFFT_COMPLEX *eigenfunction, const int m)
{
    if (test_cg==1) TITLE("Diago_CG_GPU","orthogonal_gradient");
    //timer::tick("Diago_CG_GPU","orth_grad");

    // pw.s_1psi(dim, g, sg); 
    cudaMemcpy(sg, g, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);// TODO

    int inc=1;

    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;
    // ONE ZERO cufftcomplex?
    // cublasZgemv(handle, trans1, dim, m, ONE, eigenfunction, dmx, sg, inc, ZERO, lagrange, inc);
    CUFFT_COMPLEX ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasZgemv(handle, trans1, dim, m, &ONE, eigenfunction, dmx, sg, inc, &ZERO, lagrange, inc);
    /*for (int i=0; i<m; i++)
    {
        lagrange[i] = ZERO;
        for (int j=0; j<dim; j++)
        {
            lagrange[i] += conj( eigenfunction(i,j) ) * sg[j];
        }
    }*/

    // Parallel_Reduce::reduce_complex_double_pool(lagrange, m); // todo

    // (3) orthogonal |g> and |Sg> to all states (0~m-1)

    // char trans2='N';
    cublasOperation_t trans2 = CUBLAS_OP_N;
    // cublasZgemv(handle, trans2, dim, m, NEG_ONE, eigenfunction, dmx, lagrange, inc, ONE, g, inc);
    // cublasZgemv(handle, trans2, dim, m, NEG_ONE, eigenfunction, dmx, lagrange, inc, ONE, sg, inc);

    cublasZgemv(handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, g, inc);
    cublasZgemv(handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, sg, inc);
    // TODO: 这里调了两次乘法，但事实上只需要调一次。。具体性能有待测试
    // zgemv_(&trans2,&dim,&m,&NEG_ONE,eigenfunction.c,&dmx,lagrange,&inc,&ONE,g,&inc);
    // zgemv_(&trans2,&dim,&m,&NEG_ONE,eigenfunction.c,&dmx,lagrange,&inc,&ONE,sg,&inc);
    
    /*for (int i=0; i<m; i++)
    {
        for (int j=0; j<dim; j++)
        {
            const complex<double> oo = lagrange[i] * eigenfunction(i, j);
            g[j] -= oo;
            sg[j] -= oo;
        }
    }*/

    //timer::tick("Diago_CG_GPU","orth_grad");
    cublasDestroy(handle);
    return;
}

void Diago_CG_GPU::calculate_gamma_cg(
    const int iter,
    const int dim,
    const double *precondition,
    const CUFFT_COMPLEX *g,
    const CUFFT_COMPLEX *sg,
    CUFFT_COMPLEX *psg,
    CUFFT_COMPLEX *cg,
    double &gg_last,
    const double &cg_norm,
    const double &theta,
    const CUFFT_COMPLEX *psi_m)
{
    if (test_cg==1) TITLE("Diago_CG_GPU","calculate_gamma_cg");
    //timer::tick("Diago_CG_GPU","gamma_cg");
    double gg_inter;
    if (iter>0)
    {
        // (1) Update gg_inter!
        // gg_inter = <g|psg>
        // Attention : the 'g' in psg is getted last time
        gg_inter = this->ddot_real( dim, g, psg );// b means before
    }

    // (2) Update for psg!
    // two usage:
    // firstly, for now, calculate: gg_now
    // secondly, prepare for the next iteration: gg_inter
    // |psg> = P | Sg >
    // for (int i=0; i<dim; i++)
    // {
    //     psg[i] = precondition[i] * sg[i];
    // }

    int thread = 512;
    int block = dim / thread + 1;
    kernel_precondition_inverse<<<block, thread>>>(psg, sg, dim, precondition);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

    // (3) Update gg_now!
    // gg_now = < g|P|sg > = < g|psg >
    const double gg_now = this->ddot_real( dim, g, psg );

    if (iter==0)
    {
        // (40) gg_last first value : equal gg_now
        gg_last = gg_now;
        // (50) cg direction first value : |g>
        // |cg> = |g>

        // for (int i=0; i<dim; i++)
        // {
        //     cg[i] = g[i];
        // }
        cudaMemcpy(cg, g, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
    }
    else
    {
        // (4) Update gamma !
        assert( gg_last != 0.0 );
        const double gamma = (gg_now - gg_inter) / gg_last;

        // (5) Update gg_last !
        gg_last = gg_now;

        // (6) Update cg direction !(need gamma and |go> ):
        // for (int i=0; i<dim; i++)
        // {
        //     cg[i] = gamma * cg[i] + g[i];
        // }

        kernel_get_gammacg<<<block, thread>>>(dim, cg, g, gamma);

        const double norma = gamma * cg_norm * sin(theta);
        // for (int i = 0;i < dim;i++)
        // {
        //     cg[i] -= norma * psi_m[i];
        // }

        kernel_get_normacg<<<block, thread>>>(dim, cg, psi_m, norma);
    }
    //timer::tick("Diago_CG_GPU","gamma_cg");
    return;
}


bool Diago_CG_GPU::update_psi(
    const int dim,
    double &cg_norm,
    double &theta,
    CUFFT_COMPLEX *hcg,
    const CUFFT_COMPLEX *cg,
    CUFFT_COMPLEX *scg,
    CUFFT_COMPLEX *psi_m ,
    double &eigenvalue,
    const double &threshold,
    CUFFT_COMPLEX *hpsi,
    CUFFT_COMPLEX *sphi)
{
    if (test_cg==1) TITLE("Diago_CG_GPU","update_psi");
    //timer::tick("Diago_CG_GPU","update");
    int thread = 512;
    int block = dim / 512 + 1;
    // pw.h_1psi(dim, cg, hcg, scg); // TODO
    // to cpu

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

    complex<double> *cg_cpu = new complex<double>[dim];
    complex<double> *hcg_cpu = new complex<double>[dim];
    complex<double> *scg_cpu = new complex<double>[dim];
    cudaMemcpy(cg_cpu, cg, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    cudaMemcpy(hcg_cpu, hcg, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    cudaMemcpy(scg_cpu, scg, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    
    GlobalC::hm.hpw.h_1psi(dim, cg_cpu, hcg_cpu, scg_cpu);

    err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

    // cudaMemcpy(cg, cg_cpu, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);
    cudaMemcpy(hcg, hcg_cpu, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);
    cudaMemcpy(scg, scg_cpu, dim * sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);

    delete [] cg_cpu;
    delete [] hcg_cpu;
    delete [] scg_cpu;
    
    // hpsi end
    err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

    cg_norm = sqrt( this->ddot_real(dim, cg, scg) );

    err = cudaGetLastError();
    if (err != cudaSuccess)
        cout << "Cuda error: "<< cudaGetErrorString(err) <<" in "<< __LINE__ << endl;

    if (cg_norm < 1.0e-10 ) return 1;

    const double a0 = this->ddot_real(dim, psi_m, hcg) * 2.0 / cg_norm;
    const double b0 = this->ddot_real(dim, cg, hcg) / ( cg_norm * cg_norm ) ;

    const double e0 = eigenvalue;

    theta = atan( a0/ (e0-b0) )/2.0;

    const double new_e = (e0 - b0) * cos(2.0*theta) + a0 * sin(2.0*theta);

    const double e1 = ( e0 + b0 + new_e ) /2.0;
    const double e2 = ( e0 + b0 - new_e ) /2.0;
    if (e1>e2)
    {
        theta +=  PI_HALF;
    }

    eigenvalue = min( e1, e2 );
    cout<< "====== Get E:"<<eigenvalue <<endl;

    const double cost = cos(theta);
    const double sint_norm = sin(theta)/cg_norm;

//	cout << "\n cg_norm = " << this->ddot(dim, cg, cg);
//	cout << "\n cg_norm_fac = "<< cg_norm * cg_norm;
//	cout << "\n overlap = "  << this->ddot(dim, psi_m, psi_m);

    // for (int i=0; i<dim; i++)
    // {
    //     psi_m[i] = psi_m[i] * cost + sint_norm * cg[i];
    // }

    kernel_multi_add<<<block, thread>>>(psi_m, psi_m, cost, cg, sint_norm, dim);

//	cout << "\n overlap2 = "  << this->ddot(dim, psi_m, psi_m);

    if ( abs(eigenvalue-e0)< threshold)
    {
        //timer::tick("Diago_CG_GPU","update");
        return 1;
    }
    else
    {
        // for (int i=0; i<dim; i++)
        // {
        //     sphi[i] = sphi[i] * cost + sint_norm * scg[i];
        //     hpsi[i] = hpsi[i] * cost + sint_norm * hcg[i];
        // }
        kernel_multi_add<<<block, thread>>>(sphi, sphi, cost, scg, sint_norm, dim);
        kernel_multi_add<<<block, thread>>>(hpsi, hpsi, cost, hcg, sint_norm, dim);
        //timer::tick("Diago_CG_GPU","update");
        return 0;
    }
}

// void Diago_CG_GPU::schmit_orth
// (
//     const int& dim,
//     const int& dmx,
//     const int& m,     //end
//     const CUFFT_COMPLEX *psi, // matrix
//     CUFFT_COMPLEX *sphi,
//     CUFFT_COMPLEX *psi_m
// )
// {
// //	TITLE("Diago_CG_GPU","schmit_orth");
//     //timer::tick("Diago_CG_GPU","schmit_orth");
//     // orthogonalize starting eigenfunction to those already calculated
//     // psi_m orthogonalize to psi(start) ~ psi(m-1)
//     // Attention, the orthogonalize here read as
//     // psi(m) -> psi(m) - \sum_{i < m} < psi(i) | S | psi(m) > psi(i)
//     // so the orthogonalize is performed about S.
    
//     // assert( m >= 0 );
//     cout<<"orth, dim="<<dim<<endl;
//     // assert( psi.nr >= m ); // todo
//     // CUFFT_COMPLEX *lagrange = new CUFFT_COMPLEX[ m+1 ];
//     CUFFT_COMPLEX *lagrange;
//     cudaMalloc((void**)&lagrange, (m+1)*sizeof(CUFFT_COMPLEX));

//     int inc=1;
//     int mp1 = m+1;
    
//     cublasHandle_t handle;
//     cublasCreate(&handle);
//     cublasOperation_t trans1 = CUBLAS_OP_C; 

//     CUFFT_COMPLEX ONE, ZERO, NEG_ONE;
//     ONE.y = ZERO.x = ZERO.y = 0.0;
//     ONE.x = 1.0;
//     NEG_ONE.x = -1.0;
//     cublasZgemv(handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc);

//     // zgemv_(&trans,&dim,&mp1,&ONE,psi.c,&dmx,sphi,&inc,&ZERO,lagrange,&inc);
//     // cublasZgemv(); // 得到lagrange数组，实际上是一组内积
//     //======================================================================
//     /*for (int j = 0; j <= m; j++)
//     {
//         for (int ig=0; ig < dim; ig++)
//         {
//             lagrange[j] += conj(psi( j, ig)) * sphi[ig] ;
//         }
//     }*/

//     // be careful , here reduce m+1
//     // Parallel_Reduce::reduce_complex_double_pool( lagrange, m+1 ); // todo
//     // cudaDeviceSynchronize();

//     // double psi_norm = lagrange[m].x; // lagrange[m].x
//     // 下面实现取lagrange[m]的功能
//     double psi_norm;
//     // 可以这样替代吗
//     cudaMemcpy(&psi_norm, &(lagrange[m].x), sizeof(double), cudaMemcpyDeviceToHost);
//     // double psi_norm = 1.0;

//     // handle 重复使用？
//     cublasOperation_t trans2 = CUBLAS_OP_N;
//     cublasZgemv(handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc);
//     // zgemv_(&trans2,&dim,&m,&NEG_ONE,psi.c,&dmx,lagrange,&inc,&ONE,psi_m,&inc);

//     // cout<<"psinorm: "<<psi_norm<<endl;
//     psi_norm -= ddot_real(m, lagrange, lagrange); //next
//     // cout<<"psinorm: "<<psi_norm<<endl;

//     /*for (int j = 0; j < m; j++)
//     {
//         for (int ig =0; ig < dim; ig++)
//         {
//             psi_m[ig] -= lagrange[j] * psi(j, ig);
//         }
//         psi_norm -= ( conj(lagrange[j]) * lagrange[j] ).real();
//     }*/

//     // 输出信息
//     if ( psi_norm <= 0.0)
//     {
// 		cout << " m = " << m << endl;
// 		// for(int j=0; j<=m; ++j)
// 		// {
// 		// 	cout << "\n j = " << j << " lagrange norm = " << ( conj(lagrange[j]) * lagrange[j] ).x;
// 		// }
//         cout << " in diago_cg, psi norm = " << psi_norm << endl;
// 		cout << " If you use GNU compiler, it may due to the zdotc is unavailable." << endl;
//         // WARNING_QUIT("schmit_orth","psi_norm <= 0.0");
//     }

//     psi_norm = sqrt(psi_norm); // 

//     int thread = 512;
//     int block = dim / 512 + 1;
//     kernel_normalization<<<block, thread>>>(psi_m, dim, psi_norm);

//     pw.s_1psi(dim, psi_m, sphi);

//     cublasDestroy(handle);
//     //timer::tick("Diago_CG_GPU","schmit_orth");
//     cudaFree(lagrange);
//     return ;
// }

void Diago_CG_GPU::schmit_orth
(
    const int& dim,
    const int& dmx,
    const int& m,     //end
    const CUFFT_COMPLEX *psi, // matrix
    CUFFT_COMPLEX *sphi,
    CUFFT_COMPLEX *psi_m
)
{
    assert( m >= 0 );
    // cout<<"orth, dim="<<dim<<endl;

    CUFFT_COMPLEX *lagrange;
    cudaMalloc((void**)&lagrange, (m+1)*sizeof(CUFFT_COMPLEX));
    int inc=1;
    int mp1 = m+1;
    
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C; 

    // cout<<"Test sphi ... "<<endl;
    // complex<double> *test_sphi = new complex<double>[10];
    // cudaMemcpy(test_sphi, sphi, 10*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<10;i++)
    // {
    //     cout<<test_sphi[i].real()<<" "<<test_sphi[i].imag()<<endl;
    // }
    // delete [] test_sphi;

    CUFFT_COMPLEX ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasZgemv(handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc);

    // cout<<"Test Lagrange ... "<<endl;
    // complex<double> *test_lag = new complex<double>[m+1];
    // cudaMemcpy(test_lag, lagrange, (m+1)*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<m+1;i++)
    // {
    //     cout<<test_lag[i].real()<<" "<<test_lag[i].imag()<<endl;
    // }
    // delete [] test_lag;
  
    double psi_norm;
    cudaMemcpy(&psi_norm, &lagrange[m], sizeof(double), cudaMemcpyDeviceToHost);
    cublasOperation_t trans2 = CUBLAS_OP_N;
    cublasZgemv(handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc);

    cout<<psi_norm<<endl;

    psi_norm -= ddot_real(m, lagrange, lagrange); //next

    cout<<"Psi norm before sqrt:"<<psi_norm<<endl;

    psi_norm = sqrt(psi_norm); // 

    cout<<"Psi norm after sqrt:"<<psi_norm<<endl;

    int thread = 512;
    int block = dim / 512 + 1;
    kernel_normalization<<<block, thread>>>(psi_m, dim, psi_norm);

    // pw.s_1psi(dim, psi_m, sphi);
    // TODO
    cudaMemcpy(sphi, psi_m, dim*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);

    cublasDestroy(handle);
    //timer::tick("Diago_CG_GPU","schmit_orth");
    cudaFree(lagrange);
    return ;
}


double Diago_CG_GPU::ddot_real
(
    const int &dim,
    const CUFFT_COMPLEX* psi_L,
    const CUFFT_COMPLEX* psi_R,
    const bool reduce
)
{
    int dim2=2*dim;
    cublasHandle_t handle;
    cublasCreate(&handle);
    double result;
    cublasDdot(handle, dim2, (double*)psi_L, 1, (double*)psi_R, 1, &result);
    cublasDestroy(handle);
    return result;
}

CUFFT_COMPLEX Diago_CG_GPU::ddot
(
    const int & dim,
    const CUFFT_COMPLEX * psi_L,
    const CUFFT_COMPLEX * psi_R
)
{
    // for (int i = 0; i < dim ; i++)
    // {
    //     result += conj(psi_L[i]) *  psi_R[i] ; //核函数 or cublas
    // }
    cublasHandle_t handle;
    cublasCreate(&handle);
    CUFFT_COMPLEX result;
    cublasZdotc(handle, dim, psi_L, 1, psi_R, 1, &result);
    // Parallel_Reduce::reduce_complex_double_pool( result );
    cublasDestroy(handle);
    return result;
}  // end of ddot

// this return <psi(m)|psik>
// psi的第m行和psik的内积
CUFFT_COMPLEX Diago_CG_GPU::ddot
(
    const int & dim,
    const CUFFT_COMPLEX *psi, //complex
    const int & m,
    CUFFT_COMPLEX *psik
)
{
    // assert(dim > 0) ;
    // for (int i = 0; i < dim ; i++)
    // {
    //     result += conj(psi(m, i)) *  psik[i] ;
    // }
    //有待测试 7.7上午完成这部分测试，是否可以只传一个首地址和dim即可进行运算？
    cublasHandle_t handle;
    cublasCreate(&handle);
    CUFFT_COMPLEX result;
    cublasZdotc(handle, dim, &psi[m*dim], 1, psik, 1, &result);
    // Parallel_Reduce::reduce_complex_double_pool( result );
    cublasDestroy(handle);
    return result;
}  // end of ddot


// this return <psi_L(m) | psi_R(n)>
// psiL的第m行和psiR的第n行做内积
CUFFT_COMPLEX Diago_CG_GPU::ddot
(
    const int & dim,
    const CUFFT_COMPLEX *psi_L,
    const int & m,
    const CUFFT_COMPLEX *psi_R,
    const int & n
)
{
    // assert( (dim>0) && (dim<=psi_L.nc) && (dim<=psi_R.nc) );

    // for ( int i = 0; i < dim ; i++)
    // {
    //     result += conj( psi_L(m,i) ) * psi_R(n,i) ;
    // }
    cublasHandle_t handle;
    cublasCreate(&handle);
    CUFFT_COMPLEX result;
    cublasZdotc(handle, dim, &psi_L[m*dim], 1, &psi_R[n*dim], 1, &result);
    // Parallel_Reduce::reduce_complex_double_pool( result );
    
    cublasDestroy(handle);
    return result;
} // end of ddot
