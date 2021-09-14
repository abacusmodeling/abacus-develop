#include "diago_cg_gpu.h"
#include "cuda_runtime.h"
#include "global.h"

template<class T, class T2>
int Diago_CG_GPU<T, T2>::moved = 0;

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
Diago_CG_GPU<T, T2>::Diago_CG_GPU()
{
    test_cg=0;
    cublasCreate(&diag_handle);
    // cublasCreate(&ddot_handle);
}

template<class T, class T2>
Diago_CG_GPU<T, T2>::~Diago_CG_GPU() 
{
    cublasDestroy(diag_handle);
    // cublasDestroy(ddot_handle);
}

template<class T, class T2>
void Diago_CG_GPU<T, T2>::diag
(
    T2 *phi, // matrix nband*dim
    T *e,
    const int &dim,
    const int &dmx,
    const int &n_band,
    const T *precondition,
    const T &eps,
    const int &maxter,
    const bool &reorder,
    int &notconv,
    T &avg_iter
)
{

    // cout<<"begin diago fft dim"<<GlobalC::pw.nx<<" "<<GlobalC::pw.ny<<" "<<GlobalC::pw.nz<<endl;
    // cout << &GlobalC::pw << endl;
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_GPU","diag");
    ModuleBase::timer::tick("Diago_CG_GPU","diag");

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

    // timer::tick("Diago_CG_GPU","diag");

	// Init with ZERO ...
    // T em_host = 0;

    for (int m=0; m<n_band; m++)
    {
        if (test_cg>2) GlobalV::ofs_running << "Diagonal Band : " << m << endl;

        CHECK_CUDA(cudaMemcpy(phi_m, &phi[m*dmx], dim*sizeof(T2), cudaMemcpyDeviceToDevice));

        // CHECK_CUDA(cudaMemcpy(sphi, phi_m, dim * sizeof(T2), cudaMemcpyDeviceToDevice));
        GlobalC::hm.hpw.s_1psi_gpu(dim, phi_m, sphi);

        this->schmit_orth(dim, dmx, m, phi, sphi, phi_m);

        GlobalC::hm.hpw.h_1psi_gpu(dim, phi_m, hphi, sphi);

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
            this->calculate_gradient( precondition, dim, hphi, sphi, g, pphi );
            this->orthogonal_gradient( dim, dmx, g, scg, lagrange, phi, m );
            this->calculate_gamma_cg( iter, dim, precondition, g, scg,
			    g0, cg, gg_last, cg_norm, theta, phi_m);// scg used as sg
            converged = this->update_psi( dim, cg_norm, theta, pphi, cg, scg, phi_m ,
			    em_host, eps, hphi, sphi); // pphi is used as hcg
            cudaMemcpy(&e[m], &em_host, sizeof(T), cudaMemcpyHostToDevice);
            if ( converged ) break;
        }//end iter

        CHECK_CUDA(cudaMemcpy(&phi[m*dmx], phi_m, dim*sizeof(T2), cudaMemcpyDeviceToDevice));

        if (!converged)
        {
            ++notconv;
        }

        avg_iter += static_cast<T>(iter) + 1.00;

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

    // timer::tick("Diago_CG_GPU","diag");
    CHECK_CUDA(cudaFree(lagrange));
    CHECK_CUDA(cudaFree(pphi));
    CHECK_CUDA(cudaFree(g0));
    CHECK_CUDA(cudaFree(cg));
    CHECK_CUDA(cudaFree(g));
    CHECK_CUDA(cudaFree(hphi));
    CHECK_CUDA(cudaFree(scg));
    CHECK_CUDA(cudaFree(sphi));
    CHECK_CUDA(cudaFree(phi_m));

    ModuleBase::timer::tick("Diago_CG_GPU","diag");
    return;
} // end subroutine ccgdiagg


template<class T, class T2>
void Diago_CG_GPU<T, T2>::calculate_gradient(
    const T* precondition, const int dim,
    const T2 *hpsi, const T2 *spsi,
    T2 *g, T2 *ppsi)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_GPU","calculate_gradient");
    ModuleBase::timer::tick("Diago_CG_GPU","calculate_grad");

    int thread = 512;
    int block = (dim + thread - 1) / thread;

    // kernel_precondition(data, res, size, precondition)
    // (2) PH|psi> : g[i] = hpsi[i]/precondition[i]
    kernel_precondition<<<block, thread>>>(g, hpsi, dim, precondition);
    // (3) PS|psi> : ppsi[i] = spsi[i]/precondition[i]
    kernel_precondition<<<block, thread>>>(ppsi, spsi, dim, precondition);

    // Update lambda !
    // (4) <psi|SPH|psi >
    const T eh = this->ddot_real(dim, spsi, g);
    // (5) <psi|SPS|psi >
    const T es = this->ddot_real(dim, spsi, ppsi);
    const T lambda = eh / es;

    // Update g !
    kernel_get_gredient<<<block, thread>>>(g, ppsi, dim, lambda);
    // kernel_multi_add<<<block, thread>>>(g, g, 1, ppsi, -lambda, dim);
    ModuleBase::timer::tick("Diago_CG_GPU","calculate_grad");
    return;
}


template<class T, class T2>
void Diago_CG_GPU<T, T2>::orthogonal_gradient( const int &dim, const int &dmx,
                                    T2 *g, T2 *sg, T2 *lagrange,
                                    const T2 *eigenfunction, const int m)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_GPU","orthogonal_gradient");
    ModuleBase::timer::tick("Diago_CG_GPU","orth_grad");

    GlobalC::hm.hpw.s_1psi_gpu(dim, g, sg);

    int inc=1;

    // cublasHandle_t handle;
    // cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;
    // ONE ZERO cufftcomplex?
    // cublasZgemv(handle, trans1, dim, m, ONE, eigenfunction, dmx, sg, inc, ZERO, lagrange, inc);
    T2 ONE, ZERO, NEG_ONE;
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

    ModuleBase::timer::tick("Diago_CG_GPU","orth_grad");
    // cublasDestroy(handle);
    return;
}

template<class T, class T2>
void Diago_CG_GPU<T, T2>::calculate_gamma_cg(
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
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_GPU","calculate_gamma_cg");
    ModuleBase::timer::tick("Diago_CG_GPU","gamma_cg");
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
    kernel_precondition_inverse<<<block, thread>>>(psg, sg, dim, precondition);

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

        kernel_get_gammacg<<<block, thread>>>(dim, cg, g, gamma);

        const T norma = gamma * cg_norm * sin(theta);
        // for (int i = 0;i < dim;i++)
        // {
        //     cg[i] -= norma * psi_m[i];
        // }

        kernel_get_normacg<<<block, thread>>>(dim, cg, psi_m, norma);
    }
    ModuleBase::timer::tick("Diago_CG_GPU","gamma_cg");
    return;
}


template<class T, class T2>
bool Diago_CG_GPU<T, T2>::update_psi(
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
    T2 *sphi)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG_GPU","update_psi");
    ModuleBase::timer::tick("Diago_CG_GPU","update_psi");
    int thread = 512;
    int block = (dim + thread - 1) / thread;
    // pw.h_1psi(dim, cg, hcg, scg); // TODO
    // to cpu
    GlobalC::hm.hpw.h_1psi_gpu(dim, cg, hcg, scg);
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
        theta +=  ModuleBase::PI_HALF;
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

    kernel_multi_add<<<block, thread>>>(psi_m, psi_m, cost, cg, sint_norm, dim);

//	cout << "\n overlap2 = "  << this->ddot(dim, psi_m, psi_m);

    if ( abs(eigenvalue-e0)< threshold)
    {
        ModuleBase::timer::tick("Diago_CG_GPU","update_psi");
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
        ModuleBase::timer::tick("Diago_CG_GPU","update_psi");
        return 0;
    }
}

template<class T, class T2>
void Diago_CG_GPU<T, T2>::schmit_orth
(
    const int& dim,
    const int& dmx,
    const int& m,     //end
    const T2 *psi, // matrix
    T2 *sphi,
    T2 *psi_m
)
{
    ModuleBase::timer::tick("Diago_CG_GPU","schmit_orth");
    assert( m >= 0 );
    // cout<<"orth, dim="<<dim<<endl;

    T2 *lagrange;
    CHECK_CUDA(cudaMalloc((void**)&lagrange, (m+1)*sizeof(T2)));
    int inc=1;
    int mp1 = m+1;

    // cublasHandle_t handle;
    // cublasCreate(&handle);
    cublasOperation_t trans1 = CUBLAS_OP_C;

    T2 ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    cublasZgemv(diag_handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc);

    T psi_norm;
    CHECK_CUDA(cudaMemcpy(&psi_norm, &lagrange[m], sizeof(T), cudaMemcpyDeviceToHost));
    cublasOperation_t trans2 = CUBLAS_OP_N;
    cublasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc);

    psi_norm -= ddot_real(m, lagrange, lagrange); //next
    psi_norm = sqrt(psi_norm);

    int thread = 512;
    int block = (dim + thread - 1) / thread;
    kernel_normalization<<<block, thread>>>(psi_m, dim, psi_norm);

    GlobalC::hm.hpw.s_1psi_gpu(dim, psi_m, sphi);

    // cublasDestroy(handle);
    ModuleBase::timer::tick("Diago_CG_GPU","schmit_orth");
    CHECK_CUDA(cudaFree(lagrange));
    return ;
}


template<class T, class T2>
T Diago_CG_GPU<T, T2>::ddot_real
(
    const int &dim,
    const T2* psi_L,
    const T2* psi_R,
    const bool reduce
)
{
    int dim2=2*dim;
    // cublasHandle_t handle;
    // cublasCreate(&handle);
    T result;
    // todo:
    cublasDdot(diag_handle, dim2, (T*)psi_L, 1, (T*)psi_R, 1, &result);
    // cublasDestroy(handle);
    return result;
}

template<class T, class T2>
T2 Diago_CG_GPU<T, T2>::ddot
(
    const int & dim,
    const T2 * psi_L,
    const T2 * psi_R
)
{
    // for (int i = 0; i < dim ; i++)
    // {
    //     result += conj(psi_L[i]) *  psi_R[i] ;
    // }
    // cublasHandle_t handle;
    // cublasCreate(&handle);
    T2 result;
    cublasZdotc(diag_handle, dim, psi_L, 1, psi_R, 1, &result);
    // Parallel_Reduce::reduce_complex_double_pool( result );
    // cublasDestroy(handle);
    return result;
}  // end of ddot

// this return <psi(m)|psik>
template<class T, class T2>
T2 Diago_CG_GPU<T, T2>::ddot
(
    const int & dim,
    const T2 *psi, //complex
    const int & m,
    T2 *psik
)
{
    // assert(dim > 0) ;
    // for (int i = 0; i < dim ; i++)
    // {
    //     result += conj(psi(m, i)) *  psik[i] ;
    // }
    // cublasHandle_t handle;
    // cublasCreate(&handle);
    T2 result;
    cublasZdotc(diag_handle, dim, &psi[m*dim], 1, psik, 1, &result);
    // Parallel_Reduce::reduce_complex_double_pool( result );
    // cublasDestroy(handle);
    return result;
}  // end of ddot


// this return <psi_L(m) | psi_R(n)>
template<class T, class T2>
T2 Diago_CG_GPU<T, T2>::ddot
(
    const int & dim,
    const T2 *psi_L,
    const int & m,
    const T2 *psi_R,
    const int & n
)
{
    // assert( (dim>0) && (dim<=psi_L.nc) && (dim<=psi_R.nc) );

    // for ( int i = 0; i < dim ; i++)
    // {
    //     result += conj( psi_L(m,i) ) * psi_R(n,i) ;
    // }
    // cublasHandle_t handle;
    // cublasCreate(&handle);
    T2 result;
    cublasZdotc(diag_handle, dim, &psi_L[m*dim], 1, &psi_R[n*dim], 1, &result);
    // Parallel_Reduce::reduce_complex_double_pool( result );

    // cublasDestroy(handle);
    return result;
} // end of ddot

template class Diago_CG_GPU<double, double2>;
// template class Diago_CG_GPU<float, float2>;
