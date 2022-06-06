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
    ModuleBase::timer::tick("PW_Basis_K", "real2recip");

    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.aux1[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.aux1,ft.aux1);

    this->gatherp_scatters(this->ft.aux1, this->ft.aux2);
    
    this->ft.fftzfor(ft.aux2,ft.aux1);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / double(this->nxyz) * this->ft.aux1[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.aux1[this->igl2isz_k[igl+startig]] / double(this->nxyz);
    }
    return;
    ModuleBase::timer::tick("PW_Basis_K", "real2recip");
}

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), double data
/// out: (nz, ns), complex<double> data
///
void PW_Basis_K:: real2recip(const double * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick("PW_Basis_K", "real2recip_gamma_only");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->ft.r_rspace[ir] = in[ir];
    // }
    // r2c in place
    const int npy = this->ny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy2 = ix*npy*2;
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.r_rspace[ixpy2 + ipy] = in[ixpy + ipy];
        }
    }

    this->ft.fftxyr2c(ft.r_rspace,ft.aux1);

    this->gatherp_scatters(this->ft.aux1, this->ft.aux2);
    
    this->ft.fftzfor(ft.aux2,ft.aux1);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / double(this->nxyz) * this->ft.aux1[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.aux1[this->igl2isz_k[igl+startig]] / double(this->nxyz);
    }
    ModuleBase::timer::tick("PW_Basis_K", "real2recip_gamma_only");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), complex<double>
///
void PW_Basis_K:: recip2real(const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick("PW_Basis_K", "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.aux1, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.aux1[this->igl2isz_k[igl+startig]] = in[igl];
    }
    this->ft.fftzbac(ft.aux1, ft.aux2);

    this->gathers_scatterp(this->ft.aux2,this->ft.aux1);

    this->ft.fftxybac(ft.aux1,ft.aux1);
    
    if(add)
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] += factor * this->ft.aux1[ir];
    }
    else
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.aux1[ir];
    }
    ModuleBase::timer::tick("PW_Basis_K", "recip2real");

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), double
///
void PW_Basis_K:: recip2real(const std::complex<double> * in, double * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick("PW_Basis_K", "recip2real_gamma_only");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.aux1, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.aux1[this->igl2isz_k[igl+startig]] = in[igl];
    }
   this->ft.fftzbac(ft.aux1, ft.aux2);
    
    this->gathers_scatterp(this->ft.aux2, this->ft.aux1);

    this->ft.fftxyc2r(ft.aux1,ft.r_rspace);

    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     out[ir] = this->ft.r_rspace[ir] / this->nxyz;
    // }

    // r2c in place
    const int npy = this->ny * this->nplane;
    if(add)
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy2 = ix*npy*2;
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] += factor * this->ft.r_rspace[ixpy2 + ipy];
        }
    }
    else
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy2 = ix*npy*2;
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] = this->ft.r_rspace[ixpy2 + ipy];
        }
    }
    ModuleBase::timer::tick("PW_Basis_K", "recip2real_gamma_only");
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
    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxf1[ir] = in[ir];
    }
    this->ft.fftfxyfor(ft.auxf1,ft.auxf1);

    this->gatherp_scatters(this->ft.auxf1, this->ft.auxf2);
    
    this->ft.fftfzfor(ft.auxf2,ft.auxf1);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / float(this->nxyz) * this->ft.auxf1[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxf1[this->igl2isz_k[igl+startig]] / float(this->nxyz);
    }
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), float data
/// out: (nz, ns), complex<float> data
///
void PW_Basis_K:: real2recip(const float * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    assert(this->gamma_only == true);
    const int npy = this->ny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy2 = ix*npy*2;
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.rf_rspace[ixpy2 + ipy] = in[ixpy + ipy];
        }
    }

    this->ft.fftfxyr2c(ft.rf_rspace,ft.auxf1);

    this->gatherp_scatters(this->ft.auxf1, this->ft.auxf2);
    
    this->ft.fftfzfor(ft.auxf2,ft.auxf1);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / float(this->nxyz) * this->ft.auxf1[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxf1[this->igl2isz_k[igl+startig]] / float(this->nxyz);
    }
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), complex<float>
///
void PW_Basis_K:: recip2real(const std::complex<float> * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.auxf1, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxf1[this->igl2isz_k[igl+startig]] = in[igl];
    }
    this->ft.fftfzbac(ft.auxf1, ft.auxf2);

    this->gathers_scatterp(this->ft.auxf2,this->ft.auxf1);

    this->ft.fftfxybac(ft.auxf1,ft.auxf1);
    
    if(add)
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] += factor * this->ft.auxf1[ir];
    }
    else
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.auxf1[ir];
    }

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), float
///
void PW_Basis_K:: recip2real(const std::complex<float> * in, float * out, const int ik, const bool add, const float factor)
{
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.auxf1, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxf1[this->igl2isz_k[igl+startig]] = in[igl];
    }
   this->ft.fftfzbac(ft.auxf1, ft.auxf2);
    
    this->gathers_scatterp(this->ft.auxf2, this->ft.auxf1);

    this->ft.fftfxyc2r(ft.auxf1,ft.rf_rspace);

    const int npy = this->ny * this->nplane;
    if(add)
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy2 = ix*npy*2;
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] += factor * this->ft.rf_rspace[ixpy2 + ipy];
        }
    }
    else
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy2 = ix*npy*2;
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] = this->ft.rf_rspace[ixpy2 + ipy];
        }
    }
    return;
}

#endif
}