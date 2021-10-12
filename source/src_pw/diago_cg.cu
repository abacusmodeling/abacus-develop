#include "diago_cg.cuh"
#include "cuda_runtime.h"
#include "global.h"

template<class T, class T2>
int Diago_CG_CUDA<T, T2>::moved = 0;

template<class T, class T2>
__global__ void kernel_normalization(T2 *data, int size, T norm)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        data[idx].x /= norm;
        data[idx].y /= norm;
    }
}

template<class T, class T2>
__global__ void kernel_precondition(T2 *res, const T2 *data, const int size, const T *P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        res[idx].x = data[idx].x / P[idx];
        res[idx].y = data[idx].y / P[idx];
    }
}

template<class T, class T2>
__global__ void kernel_precondition_inverse(T2 *res, const T2 *data, const int size, const T *P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        res[idx].x = data[idx].x * P[idx];
        res[idx].y = data[idx].y * P[idx];
    }
}

template<class T, class T2>
__global__ void kernel_get_gredient(T2 *g, T2 *ppsi, int size, T lambda)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        g[idx].x -= lambda * ppsi[idx].x;
        g[idx].y -= lambda * ppsi[idx].y;
    }
}

template<class T, class T2>
__global__ void kernel_get_gammacg(int size, T2 *dst, const T2 *src, T gamma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = dst[idx].x * gamma + src[idx].x;
        dst[idx].y = dst[idx].y * gamma + src[idx].y;
    }
}

template<class T, class T2>
__global__ void kernel_get_normacg(int size, T2 *dst, const T2 *src, T norma)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = dst[idx].x - norma * src[idx].x;
        dst[idx].y = dst[idx].y - norma * src[idx].y;
    }
}

template<class T, class T2>
__global__ void kernel_multi_add(T2 *dst, T2 *src1, T a1, const T2 *src2, T a2, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src1[idx].x * a1 + src2[idx].x * a2;
        dst[idx].y = src1[idx].y * a1 + src2[idx].y * a2;
    }
}

template<class T, class T2>
Diago_CG_CUDA<T, T2>::Diago_CG_CUDA()
{
    test_cg=0;
    cublasCreate(&diag_handle);
    // cublasCreate(&ddot_handle);
}

template<class T, class T2>
Diago_CG_CUDA<T, T2>::~Diago_CG_CUDA() 
{
    cublasDestroy(diag_handle);
    // cublasDestroy(ddot_handle);
}

template<class T>
void test_print(T *data, int size)
{
    T *h_data = (T*)malloc(size * sizeof(T));
    cudaMemcpy(h_data, data, size*sizeof(T), cudaMemcpyDeviceToHost);
    cout<<sizeof(h_data[0])<<endl;
    for(int i=0;i<size;i++){
        cout<<h_data[i].x<<" "<<h_data[i].y<<endl;
    }
    delete [] h_data;
}

template<class T, class T2>
void Diago_CG_CUDA<T, T2>::diag
(
    T2 *phi, // matrix nband*dim
    T *e,
    T2 *vkb_c,
    const int &dim,
    const int &dmx,
    const int &n_band,
    const T *precondition,
    const T &eps,
    const int &maxter,
    const bool &reorder,
    int &notconv,
    double &avg_iter
)
{

    // cout<<"begin diago fft dim"<<GlobalC::pw.nx<<" "<<GlobalC::pw.ny<<" "<<GlobalC::pw.nz<<endl;
    // cout << &GlobalC::pw << endl;
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_CUDA","diag");
    ModuleBase::timer::tick("Diago_CG_CUDA","diag");

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

    T2 *sphi;
    T2 *scg;
    T2 *hphi;
    T2 *g;
    T2 *cg;
    T2 *g0;
    T2 *pphi;
    T2 *lagrange;
    T2 *phi_m;

    // cout << "Hello, CG!" << endl;
    // cout << "CG Dim = " << dim << " & " << dmx << endl;

    CHECK_CUDA(cudaMalloc((void**)&sphi, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&scg, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&hphi, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&g, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&cg, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&g0, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&pphi, dim * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&lagrange, n_band * sizeof(T2)));
    CHECK_CUDA(cudaMalloc((void**)&phi_m, dim * sizeof(T2)));

    // timer::tick("Diago_CG_CUDA","diag");

	// Init with ZERO ...
    // T em_host = 0;

    for (int m=0; m<n_band; m++)
    {
        if (test_cg>2) GlobalV::ofs_running << "Diagonal Band : " << m << endl;
        // cout<<"====band====="<<m<<endl;
        cudaDeviceSynchronize();

        CHECK_CUDA(cudaMemcpy(phi_m, &phi[m*dmx], dim*sizeof(T2), cudaMemcpyDeviceToDevice));

        // CHECK_CUDA(cudaMemcpy(sphi, phi_m, dim * sizeof(T2), cudaMemcpyDeviceToDevice));
        GlobalC::hm.hpw.s_1psi_cuda(dim, phi_m, sphi);

        this->schmit_orth(dim, dmx, m, phi, sphi, phi_m);

        // cout<<"====phi_m after schmit===="<<endl;
        // test_print<T2>(phi_m, 15);

        GlobalC::hm.hpw.h_1psi_cuda(dim, phi_m, hphi, sphi, vkb_c);

        // cout<<"====hphi after hpsi===="<<endl;
        // test_print<T2>(hphi, 15);

        T em_host = 0;
        em_host = ddot_real(dim, phi_m, hphi);

        CHECK_CUDA(cudaMemcpy(&e[m], &em_host, sizeof(T), cudaMemcpyHostToDevice));

        int iter = 0;
        T gg_last = 0.0;
        T cg_norm = 0.0;
        T theta = 0.0;
        bool converged = false;
        // cg iteration

        for (iter = 0;iter < maxter; iter++)
        {
            // cout<<"******iter:"<<iter<<"******"<<endl;
            this->calculate_gradient( precondition, dim, hphi, sphi, g, pphi );

            // cout<<"====g after cal_grad===="<<endl;
            // test_print<T2>(g, 15);

            this->orthogonal_gradient( dim, dmx, g, scg, lagrange, phi, m );

            // cout<<"====lag after orth===="<<endl;
            // test_print<T2>(lagrange, 5);

            this->calculate_gamma_cg( iter, dim, precondition, g, scg,
			    g0, cg, gg_last, cg_norm, theta, phi_m);// scg used as sg
            converged = this->update_psi( dim, cg_norm, theta, pphi, cg, scg, phi_m ,
			    em_host, eps, hphi, sphi, vkb_c); // pphi is used as hcg
            
            // cout<<"====hphi after update===="<<endl;
            // test_print<T2>(hphi, 15);

            cudaMemcpy(&e[m], &em_host, sizeof(T), cudaMemcpyHostToDevice);
            if ( converged ) break;
        }//end iter

        CHECK_CUDA(cudaMemcpy(&phi[m*dmx], phi_m, dim*sizeof(T2), cudaMemcpyDeviceToDevice));

        if (!converged)
        {
            ++notconv;
        }

        // cout<<"now_iter:"<<iter<<endl;
        avg_iter += static_cast<double>(iter) + 1.00;

        if (m > 0 && reorder)
        {
			ModuleBase::GlobalFunc::NOTE("reorder bands!");
            T* e_host;
            e_host = (T*)malloc(n_band*sizeof(T));
            ModuleBase::GlobalFunc::ZEROS(e_host, n_band);
            CHECK_CUDA(cudaMemcpy(e_host, e, n_band*sizeof(T), cudaMemcpyDeviceToHost));

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
                T e0 = e_host[m];

                CHECK_CUDA(cudaMemcpy(pphi, &phi[m*dmx], dim*sizeof(T2), cudaMemcpyDeviceToDevice));

                for (int j = m;j >= i + 1;j--)
                {
                    e_host[j]=e_host[j-1];
                    CHECK_CUDA(cudaMemcpy(&phi[j*dmx], &phi[(j-1)*dmx], dim*sizeof(T2), cudaMemcpyDeviceToDevice));
                }

                e_host[i] = e0;

                CHECK_CUDA(cudaMemcpy(&phi[i*dmx], pphi, dim*sizeof(T2), cudaMemcpyDeviceToDevice));
                // this procedure should be good if only a few inversions occur,
                // extremely inefficient if eigenvectors are often in bad order
                // (but this should not happen)
            } // endif

            CHECK_CUDA(cudaMemcpy(e, e_host, n_band*sizeof(T), cudaMemcpyHostToDevice));
            delete [] e_host;
        } //end reorder

    }//end m

    avg_iter /= n_band;

    // timer::tick("Diago_CG_CUDA","diag");
    CHECK_CUDA(cudaFree(lagrange));
    CHECK_CUDA(cudaFree(pphi));
    CHECK_CUDA(cudaFree(g0));
    CHECK_CUDA(cudaFree(cg));
    CHECK_CUDA(cudaFree(g));
    CHECK_CUDA(cudaFree(hphi));
    CHECK_CUDA(cudaFree(scg));
    CHECK_CUDA(cudaFree(sphi));
    CHECK_CUDA(cudaFree(phi_m));

    ModuleBase::timer::tick("Diago_CG_CUDA","diag");
    return;
} // end subroutine ccgdiagg


template<class T, class T2>
void Diago_CG_CUDA<T, T2>::calculate_gradient(
    const T* precondition, const int dim,
    const T2 *hpsi, const T2 *spsi,
    T2 *g, T2 *ppsi)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_CUDA","calculate_gradient");
    ModuleBase::timer::tick("Diago_CG_CUDA","calculate_grad");

    int thread = 512;
    int block = (dim + thread - 1) / thread;

    // kernel_precondition(data, res, size, precondition)
    // (2) PH|psi> : g[i] = hpsi[i]/precondition[i]
    kernel_precondition<T, T2><<<block, thread>>>(g, hpsi, dim, precondition);
    // (3) PS|psi> : ppsi[i] = spsi[i]/precondition[i]
    kernel_precondition<T, T2><<<block, thread>>>(ppsi, spsi, dim, precondition);

    // Update lambda !
    // (4) <psi|SPH|psi >
    const T eh = this->ddot_real(dim, spsi, g);
    // (5) <psi|SPS|psi >
    const T es = this->ddot_real(dim, spsi, ppsi);
    const T lambda = eh / es;

    // Update g !
    kernel_get_gredient<T, T2><<<block, thread>>>(g, ppsi, dim, lambda);
    // kernel_multi_add<<<block, thread>>>(g, g, 1, ppsi, -lambda, dim);
    ModuleBase::timer::tick("Diago_CG_CUDA","calculate_grad");
    return;
}


template<class T, class T2>
void Diago_CG_CUDA<T, T2>::orthogonal_gradient( const int &dim, const int &dmx,
                                    float2 *g, float2 *sg, float2 *lagrange,
                                    const float2 *eigenfunction, const int m)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_CUDA","orthogonal_gradient");
    ModuleBase::timer::tick("Diago_CG_CUDA","orth_grad");

    GlobalC::hm.hpw.s_1psi_cuda(dim, g, sg);

    int inc=1;

    // cublasHandle_t handle;
    // cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;
    // ONE ZERO cufftcomplex?
    // cublasZgemv(handle, trans1, dim, m, ONE, eigenfunction, dmx, sg, inc, ZERO, lagrange, inc);
    float2 ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasCgemv(diag_handle, trans1, dim, m, &ONE, eigenfunction, dmx, sg, inc, &ZERO, lagrange, inc);
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
    cublasOperation_t trans2 = CUBLAS_OP_N;

    cublasCgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, g, inc);
    cublasCgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, sg, inc);

    /*for (int i=0; i<m; i++)
    {
        for (int j=0; j<dim; j++)
        {
            const complex<T> oo = lagrange[i] * eigenfunction(i, j);
            g[j] -= oo;
            sg[j] -= oo;
        }
    }*/

    ModuleBase::timer::tick("Diago_CG_CUDA","orth_grad");
    // cublasDestroy(handle);
    return;
}

template<class T, class T2>
void Diago_CG_CUDA<T, T2>::orthogonal_gradient( const int &dim, const int &dmx,
                                    double2 *g, double2 *sg, double2 *lagrange,
                                    const double2 *eigenfunction, const int m)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_CUDA","orthogonal_gradient");
    ModuleBase::timer::tick("Diago_CG_CUDA","orth_grad");

    GlobalC::hm.hpw.s_1psi_cuda(dim, g, sg);

    int inc=1;

    // cublasHandle_t handle;
    // cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;
    // ONE ZERO cufftcomplex?
    // cublasZgemv(handle, trans1, dim, m, ONE, eigenfunction, dmx, sg, inc, ZERO, lagrange, inc);
    double2 ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasZgemv(diag_handle, trans1, dim, m, &ONE, eigenfunction, dmx, sg, inc, &ZERO, lagrange, inc);
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
    cublasOperation_t trans2 = CUBLAS_OP_N;

    cublasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, g, inc);
    cublasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, sg, inc);

    /*for (int i=0; i<m; i++)
    {
        for (int j=0; j<dim; j++)
        {
            const complex<T> oo = lagrange[i] * eigenfunction(i, j);
            g[j] -= oo;
            sg[j] -= oo;
        }
    }*/

    ModuleBase::timer::tick("Diago_CG_CUDA","orth_grad");
    // cublasDestroy(handle);
    return;
}

template<class T, class T2>
void Diago_CG_CUDA<T, T2>::calculate_gamma_cg(
    const int iter,
    const int dim,
    const T *precondition,
    const T2 *g,
    const T2 *sg,
    T2 *psg,
    T2 *cg,
    T &gg_last,
    const T &cg_norm,
    const T &theta,
    const T2 *psi_m)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_CUDA","calculate_gamma_cg");
    ModuleBase::timer::tick("Diago_CG_CUDA","gamma_cg");
    T gg_inter;
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
    int block = (dim + thread - 1) / thread;
    kernel_precondition_inverse<T, T2><<<block, thread>>>(psg, sg, dim, precondition);

    // (3) Update gg_now!
    // gg_now = < g|P|sg > = < g|psg >
    const T gg_now = this->ddot_real( dim, g, psg );

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
        CHECK_CUDA(cudaMemcpy(cg, g, dim*sizeof(T2), cudaMemcpyDeviceToDevice));
    }
    else
    {
        // (4) Update gamma !
        assert( gg_last != 0.0 );
        const T gamma = (gg_now - gg_inter) / gg_last;

        // (5) Update gg_last !
        gg_last = gg_now;

        // (6) Update cg direction !(need gamma and |go> ):
        // for (int i=0; i<dim; i++)
        // {
        //     cg[i] = gamma * cg[i] + g[i];
        // }

        kernel_get_gammacg<T, T2><<<block, thread>>>(dim, cg, g, gamma);

        const T norma = gamma * cg_norm * sin(theta);
        // for (int i = 0;i < dim;i++)
        // {
        //     cg[i] -= norma * psi_m[i];
        // }

        kernel_get_normacg<T, T2><<<block, thread>>>(dim, cg, psi_m, norma);
    }
    ModuleBase::timer::tick("Diago_CG_CUDA","gamma_cg");
    return;
}


template<class T, class T2>
bool Diago_CG_CUDA<T, T2>::update_psi(
    const int dim,
    T &cg_norm,
    T &theta,
    T2 *hcg,
    const T2 *cg,
    T2 *scg,
    T2 *psi_m ,
    T &eigenvalue,
    const T &threshold,
    T2 *hpsi,
    T2 *sphi,
    T2 *vkb_c)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_CUDA","update_psi");
    ModuleBase::timer::tick("Diago_CG_CUDA","update_psi");
    int thread = 512;
    int block = (dim + thread - 1) / thread;
    // pw.h_1psi(dim, cg, hcg, scg); // TODO
    // to cpu
    GlobalC::hm.hpw.h_1psi_cuda(dim, cg, hcg, scg, vkb_c);
    // hpsi end

    cg_norm = sqrt( this->ddot_real(dim, cg, scg) );

    if (cg_norm < 1.0e-10 ) return 1;

    const T a0 = this->ddot_real(dim, psi_m, hcg) * 2.0 / cg_norm;
    const T b0 = this->ddot_real(dim, cg, hcg) / ( cg_norm * cg_norm ) ;

    const T e0 = eigenvalue;

    theta = atan( a0/ (e0-b0) )/2.0;

    const T new_e = (e0 - b0) * cos(2.0*theta) + a0 * sin(2.0*theta);

    const T e1 = ( e0 + b0 + new_e ) /2.0;
    const T e2 = ( e0 + b0 - new_e ) /2.0;
    if (e1>e2)
    {
        theta +=  (T)ModuleBase::PI_HALF;
    }

    eigenvalue = min( e1, e2 );

    const T cost = cos(theta);
    const T sint_norm = sin(theta)/cg_norm;

//	cout << "\n cg_norm = " << this->ddot(dim, cg, cg);
//	cout << "\n cg_norm_fac = "<< cg_norm * cg_norm;
//	cout << "\n overlap = "  << this->ddot(dim, psi_m, psi_m);

    // for (int i=0; i<dim; i++)
    // {
    //     psi_m[i] = psi_m[i] * cost + sint_norm * cg[i];
    // }

    kernel_multi_add<T, T2><<<block, thread>>>(psi_m, psi_m, cost, cg, sint_norm, dim);

//	cout << "\n overlap2 = "  << this->ddot(dim, psi_m, psi_m);
    // cout<<"======"<<endl;
    // cout<<abs(eigenvalue-e0)<<" "<<threshold<<endl;

    if ( abs(eigenvalue-e0)< threshold)
    {
        ModuleBase::timer::tick("Diago_CG_CUDA","update_psi");
        return 1;
    }
    else
    {
        // for (int i=0; i<dim; i++)
        // {
        //     sphi[i] = sphi[i] * cost + sint_norm * scg[i];
        //     hpsi[i] = hpsi[i] * cost + sint_norm * hcg[i];
        // }
        kernel_multi_add<T, T2><<<block, thread>>>(sphi, sphi, cost, scg, sint_norm, dim);
        kernel_multi_add<T, T2><<<block, thread>>>(hpsi, hpsi, cost, hcg, sint_norm, dim);
        ModuleBase::timer::tick("Diago_CG_CUDA","update_psi");
        return 0;
    }
}

template<class T, class T2>
void Diago_CG_CUDA<T, T2>::schmit_orth
(
    const int& dim,
    const int& dmx,
    const int& m,     //end
    const float2 *psi, // matrix
    float2 *sphi,
    float2 *psi_m
)
{
    ModuleBase::timer::tick("Diago_CG_CUDA","schmit_orth");
    assert( m >= 0 );
    // cout<<"orth, dim="<<dim<<endl;

    float2 *lagrange;
    CHECK_CUDA(cudaMalloc((void**)&lagrange, (m+1)*sizeof(float2)));
    int inc=1;
    int mp1 = m+1;

    // cublasHandle_t handle;
    // cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;

    float2 ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasCgemv(diag_handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc);

    float psi_norm;
    CHECK_CUDA(cudaMemcpy(&psi_norm, &lagrange[m], sizeof(float), cudaMemcpyDeviceToHost));
    cublasOperation_t trans2 = CUBLAS_OP_N;
    cublasCgemv(diag_handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc);

    psi_norm -= ddot_real(m, lagrange, lagrange); //next
    psi_norm = sqrt(psi_norm);

    int thread = 512;
    int block = (dim + thread - 1) / thread;
    kernel_normalization<float, float2><<<block, thread>>>(psi_m, dim, psi_norm);

    GlobalC::hm.hpw.s_1psi_cuda(dim, psi_m, sphi);

    // cublasDestroy(handle);
    ModuleBase::timer::tick("Diago_CG_CUDA","schmit_orth");
    CHECK_CUDA(cudaFree(lagrange));
    return ;
}

template<class T, class T2>
void Diago_CG_CUDA<T, T2>::schmit_orth
(
    const int& dim,
    const int& dmx,
    const int& m,     //end
    const double2 *psi, // matrix
    double2 *sphi,
    double2 *psi_m
)
{
    ModuleBase::timer::tick("Diago_CG_CUDA","schmit_orth");
    assert( m >= 0 );
    // cout<<"orth, dim="<<dim<<endl;

    double2 *lagrange;
    CHECK_CUDA(cudaMalloc((void**)&lagrange, (m+1)*sizeof(double2)));
    int inc=1;
    int mp1 = m+1;

    // cublasHandle_t handle;
    // cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;

    double2 ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasZgemv(diag_handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc);

    double psi_norm;
    CHECK_CUDA(cudaMemcpy(&psi_norm, &lagrange[m], sizeof(double), cudaMemcpyDeviceToHost));
    cublasOperation_t trans2 = CUBLAS_OP_N;
    cublasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc);

    psi_norm -= ddot_real(m, lagrange, lagrange); //next
    psi_norm = sqrt(psi_norm);

    int thread = 512;
    int block = (dim + thread - 1) / thread;
    kernel_normalization<double, double2><<<block, thread>>>(psi_m, dim, psi_norm);

    GlobalC::hm.hpw.s_1psi_cuda(dim, psi_m, sphi);

    // cublasDestroy(handle);
    ModuleBase::timer::tick("Diago_CG_CUDA","schmit_orth");
    CHECK_CUDA(cudaFree(lagrange));
    return ;
}

template<class T, class T2>
float Diago_CG_CUDA<T, T2>::ddot_real
(
    const int &dim,
    const float2* psi_L,
    const float2* psi_R,
    const bool reduce
)
{
    int dim2=2*dim;
    float result;
    cublasSdot(diag_handle, dim2, (float*)psi_L, 1, (float*)psi_R, 1, &result);
    return result;
}

template<class T, class T2>
double Diago_CG_CUDA<T, T2>::ddot_real
(
    const int &dim,
    const double2* psi_L,
    const double2* psi_R,
    const bool reduce
)
{
    int dim2=2*dim;
    double result;
    cublasDdot(diag_handle, dim2, (double*)psi_L, 1, (double*)psi_R, 1, &result);
    return result;
}


template<class T, class T2>
float2 Diago_CG_CUDA<T, T2>::ddot
(
    const int & dim,
    const float2 * psi_L,
    const float2 * psi_R
)
{
    float2 result;
    cublasCdotc(diag_handle, dim, psi_L, 1, psi_R, 1, &result);
    return result;
}  // end of ddot

template<class T, class T2>
double2 Diago_CG_CUDA<T, T2>::ddot
(
    const int & dim,
    const double2 * psi_L,
    const double2 * psi_R
)
{
    double2 result;
    cublasZdotc(diag_handle, dim, psi_L, 1, psi_R, 1, &result);
    return result;
}  // end of ddot


template<class T, class T2>
float2 Diago_CG_CUDA<T, T2>::ddot
(
    const int & dim,
    const float2 *psi, //complex
    const int & m,
    float2 *psik
)
{
    // assert(dim > 0) ;
    // for (int i = 0; i < dim ; i++)
    // {
    //     result += conj(psi(m, i)) *  psik[i] ;
    // }
    float2 result;
    cublasCdotc(diag_handle, dim, &psi[m*dim], 1, psik, 1, &result);
    return result;
}  // end of ddot

template<class T, class T2>
double2 Diago_CG_CUDA<T, T2>::ddot
(
    const int & dim,
    const double2 *psi, //complex
    const int & m,
    double2 *psik
)
{
    // assert(dim > 0) ;
    // for (int i = 0; i < dim ; i++)
    // {
    //     result += conj(psi(m, i)) *  psik[i] ;
    // }
    double2 result;
    cublasZdotc(diag_handle, dim, &psi[m*dim], 1, psik, 1, &result);
    return result;
}  // end of ddot


// this return <psi_L(m) | psi_R(n)>
template<class T, class T2>
float2 Diago_CG_CUDA<T, T2>::ddot
(
    const int & dim,
    const float2 *psi_L,
    const int & m,
    const float2 *psi_R,
    const int & n
)
{
    // assert( (dim>0) && (dim<=psi_L.nc) && (dim<=psi_R.nc) );

    // for ( int i = 0; i < dim ; i++)
    // {
    //     result += conj( psi_L(m,i) ) * psi_R(n,i) ;
    // }
    float2 result;
    cublasCdotc(diag_handle, dim, &psi_L[m*dim], 1, &psi_R[n*dim], 1, &result);
    return result;
} // end of ddot

template<class T, class T2>
double2 Diago_CG_CUDA<T, T2>::ddot
(
    const int & dim,
    const double2 *psi_L,
    const int & m,
    const double2 *psi_R,
    const int & n
)
{
    // assert( (dim>0) && (dim<=psi_L.nc) && (dim<=psi_R.nc) );

    // for ( int i = 0; i < dim ; i++)
    // {
    //     result += conj( psi_L(m,i) ) * psi_R(n,i) ;
    // }
    double2 result;
    cublasZdotc(diag_handle, dim, &psi_L[m*dim], 1, &psi_R[n*dim], 1, &result);
    return result;
} // end of ddot

template class Diago_CG_CUDA<double, double2>;
template class Diago_CG_CUDA<float, float2>;
