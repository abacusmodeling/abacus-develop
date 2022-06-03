#include "use_fft.h"
#include "global.h"
#include "../module_base/memory.h"


Use_FFT::Use_FFT()
{
	porter = new std::complex<double>[1]; //qianrui fix a bug
}

Use_FFT::~Use_FFT()
{
	delete[] porter;
#ifdef __CUDA
    cudaFree(d_porter);
    cufftDestroy(fft_handle);
#endif
}

void Use_FFT::allocate(void)
{
    ModuleBase::TITLE("Use_FFT","allocate");

    delete[] porter;
    porter = new std::complex<double>[GlobalC::rhopw->nrxx];
#ifdef __CUDA
    cufftPlan3d(&fft_handle, GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz, CUFFT_Z2Z);
    cudaMalloc((void**)&d_porter, GlobalC::rhopw->nrxx * sizeof(double2));
#endif
    ModuleBase::Memory::record("Use_FFT","porter",GlobalC::rhopw->nrxx,"complexmatrix");
    return;
}

void Use_FFT::ToRealSpace(const int &is, ModuleBase::ComplexMatrix &vg, double *vr, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->recip2real(&vg(is,0), vr);
    return;
}

// void Use_FFT::ToRealSpace_psi(const int &ik, const std::complex<double> *psig, std::complex<double> *psir)
// {
// 	ModuleBase::GlobalFunc::ZEROS(psir, GlobalC::rhopw->nrxx);
// 	for(int ig=0; ig<GlobalC::wf.npw; ig++)
// 	{
// 		psir[ GlobalC::pw.ig2fftc[ GlobalC::wf.igk(ik,ig) ] ] = psig[ig];
// 	}
// 	// 1: to real space.
// 	GlobalC::pw.FFT_wfc.FFT3D(psir, 1);

// 	return;
// }


// void Use_FFT::ToRealSpace_psi(const int &ik, const int &ib, const ModuleBase::ComplexMatrix &evc, std::complex<double> *psir)
// {
// 	// (1) set value
//     ModuleBase::GlobalFunc::ZEROS( psir, GlobalC::rhopw->nrxx );
//     for (int ig=0; ig<GlobalC::wf.npw; ig++)
//     {
//         psir[ GlobalC::pw.ig2fftc[ GlobalC::wf.igk(ik,ig) ] ] = evc(ib, ig);
//     }

// 	// (2) fft and get value
//     GlobalC::pw.FFT_wfc.FFT3D(psir, 1);
//     return;
// }



void Use_FFT::ToRealSpace(const int &is, ModuleBase::ComplexMatrix &vg, ModuleBase::matrix &vr, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->recip2real(&vg(is,0), &vr(is,0));
    return;
}


// Fourer transform of vg,
// then put vg into vr.
void Use_FFT::ToRealSpace(std::complex<double> *vg, double *vr, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->recip2real(vg, vr);
    return;
}

void Use_FFT::ToReciSpace(double* vr, std::complex<double> *vg, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->real2recip(vr, vg);
	return;
}
