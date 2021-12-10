#include "fft.h"
#include <complex>
#include "pw_basis.h"
#include <cassert>
#include "../module_base/global_function.h"
#include "../module_base/timer.h"

namespace ModulePW
{

//
//transform real space to reciprocal space
//in: (nplane,nx,ny)
//out: (nz, ns)
//
void PW_Basis:: real2recip(std::complex<double> * in, std::complex<double> * out)
{
    ModuleBase::timer::tick("PW_Basis", "real2recip");

    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.aux1[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.aux1,ft.aux1);

    this->gatherp_scatters(this->ft.aux1, this->ft.aux2);
    
    this->ft.fftzfor(ft.aux2,ft.aux1);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.aux1[this->ig2isz[ig]];
    }
    return;
    ModuleBase::timer::tick("PW_Basis", "real2recip");
}

//
//transform real space to reciprocal space
//in: (nplane,nx,ny)
//out: (nz, ns)
//
void PW_Basis:: real2recip(double * in, std::complex<double> * out)
{
    ModuleBase::timer::tick("PW_Basis", "real2recip_gamma_only");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->ft.r_rspace[ir] = in[ir];
    // }
    // r2c in place
    int npy = this->bigny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.r_rspace[ix*npy*2 + ipy] = in[ix*npy + ipy];
        }
    }

    this->ft.fftxyr2c(ft.r_rspace,ft.aux1);

    this->gatherp_scatters(this->ft.aux1, this->ft.aux2);
    
    this->ft.fftzfor(ft.aux2,ft.aux1);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.aux1[this->ig2isz[ig]];
    }
    ModuleBase::timer::tick("PW_Basis", "real2recip_gamma_only");
    return;
}

//
//transform real space to reciprocal space
//in: (nz,ns)
//out: (nplane, nx, ny)
//
void PW_Basis:: recip2real(std::complex<double> * in, std::complex<double> * out)
{
    ModuleBase::timer::tick("PW_Basis", "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.aux1, this->nst * this->nz);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.aux1[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftzbac(ft.aux1, ft.aux2);

    this->gathers_scatterp(this->ft.aux2,this->ft.aux1);

    this->ft.fftxybac(ft.aux1,ft.aux1);
    
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.aux1[ir] / double(this->bignxyz);
    }
    ModuleBase::timer::tick("PW_Basis", "recip2real");

    return;
}

//
//transform real space to reciprocal space
//in: (nz,ns)
//out: (nplane, nx, ny)
//
void PW_Basis:: recip2real(std::complex<double> * in, double * out)
{
    ModuleBase::timer::tick("PW_Basis", "recip2real_gamma_only");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.aux1, this->nst * this->nz);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.aux1[this->ig2isz[ig]] = in[ig];
    }
   this->ft.fftzbac(ft.aux1, ft.aux2);
    
    this->gathers_scatterp(this->ft.aux2, this->ft.aux1);

    this->ft.fftxyc2r(ft.aux1,ft.r_rspace);

    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     out[ir] = this->ft.r_rspace[ir] / this->bignxyz;
    // }

    // r2c in place
    int npy = this->bigny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ix*npy + ipy] = this->ft.r_rspace[ix*npy*2 + ipy] / double(this->bignxyz);
        }
    }
    ModuleBase::timer::tick("PW_Basis", "recip2real_gamma_only");
    return;
}

}