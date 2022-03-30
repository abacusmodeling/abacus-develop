#include "hamiltpw.h"

namespace ModuleHamilt
{

void HamiltPW::ch_mock()
{
    return;
}

void HamiltPW::hk_mock()
{
    return;
}

void HamiltPW::hpsi_mock(const ModulePsi::Psi<std::complex<double>>& psi, ModulePsi::Psi<std::complex<double>>& hpsi) const
{
    for(size_t iband=0; iband < hpsi.get_nbands(); ++iband )
    {
        for(size_t ibasis=0; ibasis < hpsi.get_nbasis(); ++ibasis)
        {
            hpsi(iband, ibasis) = psi(iband, ibasis);
        }
    }
}

}