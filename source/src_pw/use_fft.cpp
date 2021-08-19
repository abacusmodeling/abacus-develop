#include "use_fft.h"
#include "global.h"


Use_FFT::Use_FFT()
{
	porter = new complex<double>[GlobalC::pw.nrxx];	
}

Use_FFT::~Use_FFT()
{
	delete[] porter;	
}

void Use_FFT::allocate(void)
{
    // cout<<"before title!, pw.nrxx = "<<GlobalC::pw.nrxx<<endl;
    TITLE("Use_FFT","allocate");

    delete[] porter;
    // cout<<"before new"<<endl;
    // cout<<"pw.nrxx = "<<GlobalC::pw.nrxx<<endl;
    porter = new complex<double>[GlobalC::pw.nrxx];
    // cout<<"after new"<<endl;
    Memory::record("Use_FFT","porter",GlobalC::pw.nrxx,"complexmatrix");

    return;
}

void Use_FFT::RoundTrip(
    const complex<double> *psi,
    const double *vr,
    const int *fft_index,
    complex<double> *psic)
{
	// (1) set value
    for (int ig=0; ig< GlobalC::wf.npw; ig++)
    {
        psic[ fft_index[ig]  ] = psi[ig];
    }

	// (2) fft to real space and doing things.
    GlobalC::pw.FFT_wfc.FFT3D( psic, 1);
    for (int ir=0; ir< GlobalC::pw.nrxx; ir++)
    {
        psic[ir] *= vr[ir];
    }

	// (3) fft back to G space.
    GlobalC::pw.FFT_wfc.FFT3D( psic, -1);
    return;
}

// void Use_FFT::RoundTrip_GPU(
//     const CUFFT_COMPLEX *psi,
//     const double *vr,
//     const int *fft_index,
//     CUFFT_COMPLEX *psic)
// {
//     cout<<"rounftrip on GPU!"<<endl;
//     // (1) set value
//     int thread = 512;
//     int block = GlobalC::wf.npw / thread + 1;
//     // kernel_set<<<block, thread>>>(GlobalC::wf.npw, psic, psi, fft_index);
//     // for(int ig=0;ig<wf.npw;ig++)
//     // {
//     //     psic[fft_index[ig]] = psi[ig];
//     // }

//     // (2) fft to real space and do things.
//     cufftHandle cufftplan_gpu;
//     cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
//     cufftExecZ2Z(cufftplan_gpu, psic, psic, CUFFT_FORWARD);
    
//     int block2 = GlobalC::pw.nrxx / thread + 1;
//     // kernel_roundtrip<<<block2, thread>>>(GlobalC::pw.nrxx, psic, vr);

//     // (3) fft back to G space
//     cufftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
//     cufftExecZ2Z(cufftplan_gpu, psic, psic, CUFFT_INVERSE);

//     cufftDestroy(cufftplan_gpu);

//     cout<<"rounftrip end"<<endl;
//     return;
// }

void Use_FFT::ToRealSpace(const int &is, const ComplexMatrix &vg, double *vr)
{
	// (1) set value
    ZEROS( porter, GlobalC::pw.nrxx );
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        porter[ GlobalC::pw.ig2fftc[ig] ] = vg(is, ig);
    }

	// (2) fft and get value
    GlobalC::pw.FFT_chg.FFT3D(porter, 1);
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        vr[ir] = porter[ir].real();
    }
    return;
}

void Use_FFT::ToRealSpace_psi(const int &ik, const complex<double> *psig, complex<double> *psir)
{
	ZEROS(psir, GlobalC::pw.nrxx);
	for(int ig=0; ig<GlobalC::wf.npw; ig++)
	{
		psir[ GlobalC::pw.ig2fftc[ GlobalC::wf.igk(ik,ig) ] ] = psig[ig];
	}
	// 1: to real space.
	GlobalC::pw.FFT_wfc.FFT3D(psir, 1);

	return;
}


void Use_FFT::ToRealSpace_psi(const int &ik, const int &ib, const ComplexMatrix &evc, complex<double> *psir)
{
	// (1) set value
    ZEROS( psir, GlobalC::pw.nrxx );
    for (int ig=0; ig<GlobalC::wf.npw; ig++)
    {
        psir[ GlobalC::pw.ig2fftc[ GlobalC::wf.igk(ik,ig) ] ] = evc(ib, ig);
    }

	// (2) fft and get value
    GlobalC::pw.FFT_wfc.FFT3D(psir, 1);
    return;
}



void Use_FFT::ToRealSpace(const int &is, const ComplexMatrix &vg, matrix &vr)
{
	// (1) set value
    ZEROS( porter, GlobalC::pw.nrxx);
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        porter [GlobalC::pw.ig2fftc[ig]] = vg(is,ig);
    }

	// (2) fft and get value
    GlobalC::pw.FFT_chg.FFT3D(porter, 1);
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        vr(is,ir) = porter[ir].real();
    }
    return;
}


// Fourer transform of vg, 
// then put vg into vr.
void Use_FFT::ToRealSpace(const complex<double> *vg, double *vr)
{
    ZEROS( porter, GlobalC::pw.nrxx);
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        porter[GlobalC::pw.ig2fftc[ig]] = vg[ig];
    }
    GlobalC::pw.FFT_chg.FFT3D(porter, 1);
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        vr[ir] = porter[ir].real();
    }
    return;
}

void Use_FFT::ToReciSpace(const double* vr, complex<double> *vg)
{
	ZEROS( porter, GlobalC::pw.nrxx);
	for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		porter[ir] = complex<double>(vr[ir], 0.0);
	}
	GlobalC::pw.FFT_chg.FFT3D( porter, -1);
	for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
	{
		vg[ig] = porter[ GlobalC::pw.ig2fftc[ig] ];
	}
	return;
}

