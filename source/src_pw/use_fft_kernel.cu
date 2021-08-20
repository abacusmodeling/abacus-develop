#include "use_fft_kernel.h"
#include "global.h"
#include "cufft.h"

__global__ void kernel_set(int size, CUFFT_COMPLEX *dst, const CUFFT_COMPLEX *src, const int *index_list)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int p = index_list[idx];
    if(idx < size)
    {
        dst[p].x = src[idx].x;
        dst[p].y = src[idx].y;
    }
}

__global__ void kernel_roundtrip(int size, CUFFT_COMPLEX *dst, const double *src)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x *= src[idx];
        dst[idx].y *= src[idx];
    }
}

__global__ void kernel_normalization(int size, CUFFT_COMPLEX *data, double norm)
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

void Reorder_psi_plus(CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src)
{
    int ii = 0;
    int size_z = GlobalC::pw.FFT_wfc.npps[0];
    int thread = 512;
    int block = size_z / 512 + 1;
    for(int is=0; is<GlobalC::pw.FFT_wfc.nst; is++)
    {
        int ir = GlobalC::pw.FFT_wfc.ismap[is];
        kernel_reorder<<<block, thread>>>(&dst[ir*size_z], &src[ii*size_z], size_z);
        // cudaMemcpy(&dst[ir*size_z], &src[ii*size_z], size_z*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
        ii++;
    }
}

void Reorder_psi_minus(CUFFT_COMPLEX *dst, CUFFT_COMPLEX *src)
{
    int ii = 0;
    int size_z = GlobalC::pw.FFT_wfc.npps[0];
    int thread = 512;
    int block = size_z / 512 + 1;
    for(int j=0; j<GlobalC::pw.FFT_wfc.nst; j++)
    {
        int ir = GlobalC::pw.FFT_wfc.ismap[j];
        kernel_reorder<<<block, thread>>>(&dst[ii*size_z], &src[ir*size_z], size_z);
        // cudaMemcpy(&dst[ii*size_z], &src[ir*size_z], size_z*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToDevice);
        ii++;
    }
}

void RoundTrip(const CUFFT_COMPLEX *psi, const double *vr, const int *fft_index, CUFFT_COMPLEX *psic)
{
    // (1) set value
    int thread = 512;
    int block = GlobalC::wf.npw / thread + 1;
    int block2 = GlobalC::pw.nrxx / thread + 1;
    kernel_set<<<block, thread>>>(GlobalC::wf.npw, psic, psi, fft_index);

    // todo : reorder psi
    // int *d_ismap;
    // cudaMalloc((void**)&d_ismap, GlobalC::pw.FFT_wfc.nst*sizeof(int));
    // cudaMemcpy(d_ismap, GlobalC::pw.FFT_wfc.ismap, GlobalC::pw.FFT_wfc.nst*sizeof(int), cudaMemcpyHostToDevice);

    CUFFT_COMPLEX *ordered_psi;
    cudaMalloc((void**)&ordered_psi, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX));
    cudaMemset(ordered_psi, 0, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX));

    Reorder_psi_plus(ordered_psi, psic);

    cufftHandle cufftplan_gpu;
    cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cufftExecZ2Z(cufftplan_gpu, ordered_psi, ordered_psi, CUFFT_INVERSE);
    cufftDestroy(cufftplan_gpu);

    // int block3 = GlobalC::pw.nrxx / thread + 1;
    // kernel_normalization<<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    kernel_roundtrip<<<block2, thread>>>(GlobalC::pw.nrxx, ordered_psi, vr);

    cufftHandle cufftplan_gpu2;
    cufftPlan3d(&cufftplan_gpu2, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cufftExecZ2Z(cufftplan_gpu2, ordered_psi, ordered_psi, CUFFT_FORWARD);
    cufftDestroy(cufftplan_gpu2);

    Reorder_psi_minus(psic, ordered_psi);

    int block3 = GlobalC::pw.nrxx / thread + 1;
    kernel_normalization<<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    cudaFree(ordered_psi);
/*
    // todo: delete!
    // cout<<"size of ismap"<<GlobalC::pw.FFT_wfc.nst<<endl;
    // cout<<"nppz:"<<GlobalC::pw.FFT_wfc.npps[0]<<endl;



    // cout<<"before fft1"<<endl;
    // complex<double> *psic_inside = new complex<double>[15];
    // cudaMemcpy(psic_inside, &psic[6000], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<15;i++)
    // {
    //     cout<<psic_inside[i].real()<<" "<<psic_inside[i].imag()<<endl;
    // }
    // // delete [] psic_inside;
    // cout<<"========"<<endl;


    // for(int ig=0;ig<wf.npw;ig++)
    // {
    //     psic[fft_index[ig]] = psi[ig];
    // }


    // cufftHandle cufftplan_gpu;
    // cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    // cufftExecZ2Z(cufftplan_gpu, psic, psic, CUFFT_INVERSE);
    // cufftDestroy(cufftplan_gpu);

    // complex<double> *psic_cpu = new complex<double>[GlobalC::pw.nrxx];
    // cudaMemcpy(psic_cpu, psic, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // GlobalC::pw.FFT_wfc.FFT3D( psic_cpu, 1);

    // cudaMemcpy(psic, psic_cpu, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);

    // cout<<"after fft1"<<endl;
    // psic_inside = new complex<double>[15];
    // cudaMemcpy(psic_inside, &psic[6000], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<15;i++)
    // {
    //     cout<<psic_inside[i].real()<<" "<<psic_inside[i].imag()<<endl;
    // }
    // // delete [] psic_inside;
    // cout<<"========"<<endl;

    // double *vr_cpu = new double[15];
    // cudaMemcpy(vr_cpu, &vr[6000], sizeof(double)*15, cudaMemcpyDeviceToHost);
    // cout<<"vr ERROR:"<<endl;
    // for(int i=0;i<15;i++){
    //     cout<<vr_cpu[i]<<endl;
    // }
    // cout<<"==========="<<endl;

    // kernel_roundtrip<<<block2, thread>>>(GlobalC::pw.nrxx, psic, vr);

    // cout<<"before fft2"<<endl;
    // psic_inside = new complex<double>[15];
    // cudaMemcpy(psic_inside, &psic[6000], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<15;i++)
    // {
    //     cout<<psic_inside[i].real()<<" "<<psic_inside[i].imag()<<endl;
    // }
    // delete [] psic_inside;
    // cout<<"========"<<endl;

    // (3) fft back to G space
    // cufftHandle cufftplan_gpu2;
    // cufftPlan3d(&cufftplan_gpu2, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    // cufftExecZ2Z(cufftplan_gpu2, psic, psic, CUFFT_FORWARD);

    // cufftDestroy(cufftplan_gpu2);

    // cout<<"before normalization"<<endl;
    // complex<double> *tmp2 = new complex<double>[15];
    // cudaMemcpy(tmp2, &psic[6000], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<15;i++)
    // {
    //     cout<<tmp2[i].real()<<" "<<tmp2[i].imag()<<endl;
    // }
    // delete [] tmp2;

    // int block3 = GlobalC::pw.nrxx / thread + 1;
    // kernel_normalization<<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    // complex<double> *psic_cpu = new complex<double>[GlobalC::pw.nrxx];
    // psic_cpu = new complex<double>[GlobalC::pw.nrxx];
    // cudaMemcpy(psic_cpu, psic, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // GlobalC::pw.FFT_wfc.FFT3D( psic_cpu, -1);
    // cudaMemcpy(psic, psic_cpu, GlobalC::pw.nrxx*sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);

    // cout<<"after 2nd ufft SUCCESS"<<endl;
    // complex<double> *tmp1 = new complex<double>[15];
    // cudaMemcpy(tmp1, &psic[6000], 15*sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
    // for(int i=0;i<15;i++)
    // {
    //     cout<<tmp1[i].real()<<" "<<tmp1[i].imag()<<endl;
    // }
    // delete [] tmp1;

    // cout<<"rounftrip end"<<endl;

    // cout<<"fft dim: "<<GlobalC::pw.nx<<" "<<GlobalC::pw.ny<<" "<<GlobalC::pw.nz<<endl;
*/
    return;
}