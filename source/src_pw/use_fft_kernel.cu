#include "use_fft.h"
#include "global.h"
#include "cufft.h"

template<class T2>
__global__ void kernel_set(int size, T2 *dst, const T2 *src, const int *index_list)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int p = index_list[idx];
    if(idx < size)
    {
        dst[p].x = src[idx].x;
        dst[p].y = src[idx].y;
    }
}

template<class T, class T2>
__global__ void kernel_roundtrip(int size, T2 *dst, const T *src)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x *= src[idx];
        dst[idx].y *= src[idx];
    }
}

template<class T, class T2>
__global__ void kernel_normalization(int size, T2 *data, T norm)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        data[idx].x /= norm;
        data[idx].y /= norm;
    }
}

__global__ void kernel_reorder(CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x;
        dst[idx].y = src[idx].y;
    }
}

// Donghs fix 2021.9.8
// void Reorder_psi_plus(CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src)
// {
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_plus");
//     int ii = 0;
//     int size_z = GlobalC::pw.FFT_wfc.npps[0];
//     int thread = 512;
//     int block = (size_z + thread - 1) / thread;
//     for(int is=0; is<GlobalC::pw.FFT_wfc.nst; is++)
//     {
//         int ir = GlobalC::pw.FFT_wfc.ismap[is];
//         kernel_reorder<<<block, thread>>>(&dst[ir*size_z], &src[ii*size_z], size_z);
//         // CHECK_CUDA(cudaMemcpy(&dst[ir*size_z], &src[ii*size_z], size_z*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice));
//         ii++;
//     }
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_plus");
// }

// void Reorder_psi_minus(CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src)
// {
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_minus");
//     int ii = 0;
//     int size_z = GlobalC::pw.FFT_wfc.npps[0];
//     int thread = 512;
//     int block = (size_z + thread - 1) / thread;
//     for(int j=0; j<GlobalC::pw.FFT_wfc.nst; j++)
//     {
//         int ir = GlobalC::pw.FFT_wfc.ismap[j];
//         kernel_reorder<<<block, thread>>>(&dst[ii*size_z], &src[ir*size_z], size_z);
//         // CHECK_CUDA(cudaMemcpy(&dst[ii*size_z], &src[ir*size_z], size_z*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice));
//         ii++;
//     }
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_minus");
// }


void RoundTrip_kernel(const float2 *psi, const float *vr, const int *fft_index, float2 *psic)
{
    // (1) set value
    int thread = 512;
    int block = (GlobalC::wf.npw + thread - 1) / thread;
    int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_set<float2><<<block, thread>>>(GlobalC::wf.npw, psic, psi, fft_index);

    // CUFFT_COMPLEX *ordered_psi;
    // CHECK_CUDA(cudaMalloc((void**)&ordered_psi, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX)));
    // CHECK_CUDA(cudaMemset(ordered_psi, 0, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX)));

    // Reorder_psi_plus(ordered_psi, psic);

    // cufftHandle cufftplan_gpu;
    // cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cufftExecC2C(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_INVERSE);
    // cufftDestroy(cufftplan_gpu);

    // int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    // kernel_normalization<<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    kernel_roundtrip<float, float2><<<block2, thread>>>(GlobalC::pw.nrxx, psic, vr);

    // cufftHandle cufftplan_gpu2;
    // cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cufftExecC2C(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_FORWARD);
    // cufftDestroy(cufftplan_gpu);

    // Reorder_psi_minus(psic, ordered_psi);

    int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_normalization<float, float2><<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    // CHECK_CUDA(cudaFree(ordered_psi));

    return;
}

void RoundTrip_kernel(const double2 *psi, const double *vr, const int *fft_index, double2 *psic)
{
    // (1) set value
    int thread = 512;
    int block = (GlobalC::wf.npw + thread - 1) / thread;
    int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_set<double2><<<block, thread>>>(GlobalC::wf.npw, psic, psi, fft_index);

    // CUFFT_COMPLEX *ordered_psi;
    // CHECK_CUDA(cudaMalloc((void**)&ordered_psi, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX)));
    // CHECK_CUDA(cudaMemset(ordered_psi, 0, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX)));

    // Reorder_psi_plus(ordered_psi, psic);

    // cufftHandle cufftplan_gpu;
    // cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cufftExecZ2Z(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_INVERSE);
    // cufftDestroy(cufftplan_gpu);

    // int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    // kernel_normalization<<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    kernel_roundtrip<double, double2><<<block2, thread>>>(GlobalC::pw.nrxx, psic, vr);

    // cufftHandle cufftplan_gpu2;
    // cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cufftExecZ2Z(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_FORWARD);
    // cufftDestroy(cufftplan_gpu);

    // Reorder_psi_minus(psic, ordered_psi);

    int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_normalization<double, double2><<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    // CHECK_CUDA(cudaFree(ordered_psi));

    return;
}