#include "hip/hip_runtime.h"
#include "use_fft.h"
#include "global.h"
#include "hipfft.h"

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

__global__ void kernel_reorder(double2 *dst, double2 *src, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x;
        dst[idx].y = src[idx].y;
    }
}

// Donghs fix 2021.9.8
// void Reorder_psi_plus(double2 *dst, double2 *src)
// {
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_plus");
//     int ii = 0;
//     int size_z = GlobalC::pw.FFT_wfc.npps[0];
//     int thread = 512;
//     int block = (size_z + thread - 1) / thread;
//     for(int is=0; is<GlobalC::pw.FFT_wfc.nst; is++)
//     {
//         int ir = GlobalC::pw.FFT_wfc.ismap[is];
//         hipLaunchKernelGGL(kernel_reorder, dim3(block), dim3(thread), 0, 0, &dst[ir*size_z], &src[ii*size_z], size_z);
//         // CHECK_CUDA(hipMemcpy(&dst[ir*size_z], &src[ii*size_z], size_z*sizeof(double2), hipMemcpyDeviceToDevice));
//         ii++;
//     }
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_plus");
// }

// void Reorder_psi_minus(double2 *dst, double2 *src)
// {
//     ModuleBase::timer::tick("Use_FFT","reorder_psi_minus");
//     int ii = 0;
//     int size_z = GlobalC::pw.FFT_wfc.npps[0];
//     int thread = 512;
//     int block = (size_z + thread - 1) / thread;
//     for(int j=0; j<GlobalC::pw.FFT_wfc.nst; j++)
//     {
//         int ir = GlobalC::pw.FFT_wfc.ismap[j];
//         hipLaunchKernelGGL(kernel_reorder, dim3(block), dim3(thread), 0, 0, &dst[ii*size_z], &src[ir*size_z], size_z);
//         // CHECK_CUDA(hipMemcpy(&dst[ii*size_z], &src[ir*size_z], size_z*sizeof(double2), hipMemcpyDeviceToDevice));
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
    hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_set<float2>), dim3(block), dim3(thread), 0, 0, GlobalC::wf.npw, psic, psi, fft_index);

    // double2 *ordered_psi;
    // CHECK_CUDA(hipMalloc((void**)&ordered_psi, GlobalC::pw.nrxx*sizeof(double2)));
    // CHECK_CUDA(hipMemset(ordered_psi, 0, GlobalC::pw.nrxx*sizeof(double2)));

    // Reorder_psi_plus(ordered_psi, psic);

    // hipfftHandle cufftplan_gpu;
    // CHECK_CUFFT(hipfftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
    CHECK_CUFFT(hipfftExecC2C(GlobalC::UFFT.fft_handle, psic, psic, HIPFFT_BACKWARD));
    hipDeviceSynchronize();
    // CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

    // int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    // hipLaunchKernelGGL(kernel_normalization, dim3(block3), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_roundtrip<float, float2>), dim3(block2), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic, vr);

    // hipfftHandle cufftplan_gpu2;
    // CHECK_CUFFT(hipfftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
    CHECK_CUFFT(hipfftExecC2C(GlobalC::UFFT.fft_handle, psic, psic, HIPFFT_FORWARD));
    hipDeviceSynchronize();
    // CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

    // Reorder_psi_minus(psic, ordered_psi);

    int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_normalization<float, float2>), dim3(block3), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    // CHECK_CUDA(hipFree(ordered_psi));

    return;
}

void RoundTrip_kernel(const double2 *psi, const double *vr, const int *fft_index, double2 *psic)
{
    // (1) set value
    int thread = 512;
    int block = (GlobalC::wf.npw + thread - 1) / thread;
    int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_set<double2>), dim3(block), dim3(thread), 0, 0, GlobalC::wf.npw, psic, psi, fft_index);

    // cout<<"In fft"<<endl;
    // double2 *hos_psi = (double2*)malloc(10*sizeof(double2));
    // double2 *hos_psic = (double2*)malloc(10*sizeof(double2)); 
    // CHECK_CUDA(hipMemcpy(hos_psi, psi, 10*sizeof(double2), hipMemcpyDeviceToHost));
    // CHECK_CUDA(hipMemcpy(hos_psic, psic, 10*sizeof(double2), hipMemcpyDeviceToHost));
    
    // cout<<"hospsi"<<endl;
    // for(int i=0;i<10;i++){
    //     cout<<hos_psi[i].x<<" "<<hos_psi[i].y<<endl;
    // }
    // cout<<"hospsi-c"<<endl;
    // for(int i=0;i<10;i++){
    //     cout<<hos_psic[i].x<<" "<<hos_psic[i].y<<endl;
    // }
    // delete [] hos_psi;
    // delete [] hos_psic;
    
    // double2 *ordered_psi;
    // CHECK_CUDA(hipMalloc((void**)&ordered_psi, GlobalC::pw.nrxx*sizeof(double2)));
    // CHECK_CUDA(hipMemset(ordered_psi, 0, GlobalC::pw.nrxx*sizeof(double2)));

    // Reorder_psi_plus(ordered_psi, psic);

    // hipfftHandle cufftplan_gpu;
    // CHECK_CUFFT(hipfftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
    CHECK_CUFFT(hipfftExecZ2Z(GlobalC::UFFT.fft_handle, psic, psic, HIPFFT_BACKWARD));
    hipDeviceSynchronize();
    // CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

    // int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    // hipLaunchKernelGGL(kernel_normalization, dim3(block3), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    // double2 *hos_psic_aft = (double2*)malloc(10*sizeof(double2)); 
    // CHECK_CUDA(hipMemcpy(hos_psic_aft, psic, 10*sizeof(double2), hipMemcpyDeviceToHost));
    // cout<<"hospsi-c after fft"<<endl;
    // for(int i=0;i<10;i++){
    //     cout<<hos_psic_aft[i].x<<" "<<hos_psic_aft[i].y<<endl;
    // }
    // delete [] hos_psic_aft;

    hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_roundtrip<double, double2>), dim3(block2), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic, vr);

    // hipfftHandle cufftplan_gpu2;
    // CHECK_CUFFT(hipfftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
    CHECK_CUFFT(hipfftExecZ2Z(GlobalC::UFFT.fft_handle, psic, psic, HIPFFT_FORWARD));
    hipDeviceSynchronize();
    // CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

    // Reorder_psi_minus(psic, ordered_psi);

    int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_normalization<double, double2>), dim3(block3), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    // CHECK_CUDA(hipFree(ordered_psi));

    return;
}