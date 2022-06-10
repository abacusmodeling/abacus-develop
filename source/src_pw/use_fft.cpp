#include "use_fft.h"
#include "global.h"
#include "../module_base/memory.h"


Use_FFT::Use_FFT()
{
}

Use_FFT::~Use_FFT()
{
}
void Use_FFT::ToRealSpace(const int &is, const ModuleBase::ComplexMatrix &vg, double *vr, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->recip2real(&vg(is,0), vr);
    return;
}



void Use_FFT::ToRealSpace(const int &is, const ModuleBase::ComplexMatrix &vg, ModuleBase::matrix &vr, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->recip2real(&vg(is,0), &vr(is,0));
    return;
}


// Fourer transform of vg,
// then put vg into vr.
void Use_FFT::ToRealSpace(const std::complex<double> *vg, double *vr, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->recip2real(vg, vr);
    return;
}

void Use_FFT::ToReciSpace(const double* vr, std::complex<double> *vg, ModulePW::PW_Basis* rho_basis)
{
    rho_basis->real2recip(vr, vg);
	return;
}
