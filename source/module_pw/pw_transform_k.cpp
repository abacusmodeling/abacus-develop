#include "fft.h"
#include <complex>
#include "pw_basis_k.h"
#include <cassert>
#include "../module_base/global_function.h"
#include "../module_base/timer.h"
#include "pw_gatherscatter.h"

namespace ModulePW
{

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), complex<double> data
/// out: (nz, ns),  complex<double> data
///
void PW_Basis_K:: real2recip(const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");

    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxr[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.auxr,ft.auxr);

    this->gatherp_scatters(this->ft.auxr, this->ft.auxg);
    
    this->ft.fftzfor(ft.auxg,ft.auxg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / double(this->nxyz) * this->ft.auxg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxg[this->igl2isz_k[igl+startig]] / double(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), double data
/// out: (nz, ns), complex<double> data
///
void PW_Basis_K:: real2recip(const double * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->ft.r_rspace[ir] = in[ir];
    // }
    // r2c in place
    const int npy = this->ny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.r_rspace[ixpy + ipy] = in[ixpy + ipy];
        }
    }

    this->ft.fftxyr2c(ft.r_rspace,ft.auxr);

    this->gatherp_scatters(this->ft.auxr, this->ft.auxg);
    
    this->ft.fftzfor(ft.auxg,ft.auxg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / double(this->nxyz) * this->ft.auxg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxg[this->igl2isz_k[igl+startig]] / double(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), complex<double>
///
void PW_Basis_K:: recip2real(const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.auxg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxg[this->igl2isz_k[igl+startig]] = in[igl];
    }
    this->ft.fftzbac(ft.auxg, ft.auxg);

    this->gathers_scatterp(this->ft.auxg,this->ft.auxr);

    this->ft.fftxybac(ft.auxr,ft.auxr);
    
    if(add)
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] += factor * this->ft.auxr[ir];
    }
    else
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.auxr[ir];
    }
    ModuleBase::timer::tick(this->classname, "recip2real");

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), double
///
void PW_Basis_K:: recip2real(const std::complex<double> * in, double * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.auxg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxg[this->igl2isz_k[igl+startig]] = in[igl];
    }
   this->ft.fftzbac(ft.auxg, ft.auxg);
    
    this->gathers_scatterp(this->ft.auxg, this->ft.auxr);

    this->ft.fftxyc2r(ft.auxr,ft.r_rspace);

    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     out[ir] = this->ft.r_rspace[ir] / this->nxyz;
    // }

    // r2c in place
    const int npy = this->ny * this->nplane;
    if(add)
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] += factor * this->ft.r_rspace[ixpy + ipy];
        }
    }
    else
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] = this->ft.r_rspace[ixpy + ipy];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

#ifdef __MIX_PRECISION
///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), complex<float> data
/// out: (nz, ns),  complex<float> data
///
void PW_Basis_K:: real2recip(const std::complex<float> * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxfr[ir] = in[ir];
    }
    this->ft.fftfxyfor(ft.auxfr,ft.auxfr);

    this->gatherp_scatters(this->ft.auxfr, this->ft.auxfg);
    
    this->ft.fftfzfor(ft.auxfg,ft.auxfg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / float(this->nxyz) * this->ft.auxfg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxfg[this->igl2isz_k[igl+startig]] / float(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), float data
/// out: (nz, ns), complex<float> data
///
void PW_Basis_K:: real2recip(const float * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == true);
    const int npy = this->ny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.rf_rspace[ixpy + ipy] = in[ixpy + ipy];
        }
    }

    this->ft.fftfxyr2c(ft.rf_rspace,ft.auxfr);

    this->gatherp_scatters(this->ft.auxfr, this->ft.auxfg);
    
    this->ft.fftfzfor(ft.auxfg,ft.auxfg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / float(this->nxyz) * this->ft.auxfg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxfg[this->igl2isz_k[igl+startig]] / float(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), complex<float>
///
void PW_Basis_K:: recip2real(const std::complex<float> * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.auxfg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxfg[this->igl2isz_k[igl+startig]] = in[igl];
    }
    this->ft.fftfzbac(ft.auxfg, ft.auxfg);

    this->gathers_scatterp(this->ft.auxfg,this->ft.auxfr);

    this->ft.fftfxybac(ft.auxfr,ft.auxfr);
    
    if(add)
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] += factor * this->ft.auxfr[ir];
    }
    else
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.auxfr[ir];
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), float
///
void PW_Basis_K:: recip2real(const std::complex<float> * in, float * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.auxfg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxfg[this->igl2isz_k[igl+startig]] = in[igl];
    }
   this->ft.fftfzbac(ft.auxfg, ft.auxfg);
    
    this->gathers_scatterp(this->ft.auxfg, this->ft.auxfr);

    this->ft.fftfxyc2r(ft.auxfr,ft.rf_rspace);

    const int npy = this->ny * this->nplane;
    if(add)
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] += factor * this->ft.rf_rspace[ixpy + ipy];
        }
    }
    else
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] = this->ft.rf_rspace[ixpy + ipy];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

#endif
}