#include "fft.h"
#include <complex>
#include "pw_basis.h"
#include <cassert>
#include "../module_base/global_function.h"
#include "../module_base/timer.h"
#include "pw_gatherscatter.h"

namespace ModulePW {

/// 
/// transform real space to reciprocal space
/// c(k,g)=\int dr*f(r)*exp(-ig*r)
/// c(k,g)=c_k(g)*exp(ik*r)
/// c_k(g)=\int dr*f(r)*exp(-i(g+k)*r)
/// Here we calculate c(k,g)
/// in: (nplane,ny,nx), complex<double> data
/// out: (nz, ns),  complex<double> data
///
template <typename FPTYPE>
void PW_Basis:: real2recip(const std::complex<FPTYPE> * in, std::complex<FPTYPE> * out, const bool add, const FPTYPE factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");

    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.get_auxr_data<FPTYPE>()[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.get_auxr_data<FPTYPE>(),ft.get_auxr_data<FPTYPE>());

    this->gatherp_scatters(this->ft.get_auxr_data<FPTYPE>(), this->ft.get_auxg_data<FPTYPE>());
    
    this->ft.fftzfor(ft.get_auxg_data<FPTYPE>(),ft.get_auxg_data<FPTYPE>());

    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += factor / FPTYPE(this->nxyz) * this->ft.get_auxg_data<FPTYPE>()[this->ig2isz[ig]];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = this->ft.get_auxg_data<FPTYPE>()[this->ig2isz[ig]] / FPTYPE(this->nxyz);
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
}

///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), double data
/// out: (nz, ns), complex<double> data
///
template <typename FPTYPE>
void PW_Basis:: real2recip(const FPTYPE * in, std::complex<FPTYPE> * out, const bool add, const FPTYPE factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    if(this->gamma_only)
    {
        const int npy = this->ny * this->nplane;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ix = 0 ; ix < this->nx ; ++ix)
        {
            for(int ipy = 0 ; ipy < npy ; ++ipy)
            {
                this->ft.get_rspace_data<FPTYPE>()[ix*npy + ipy] = in[ix*npy + ipy];
            }
        }

        this->ft.fftxyr2c(ft.get_rspace_data<FPTYPE>(),ft.get_auxr_data<FPTYPE>());
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            this->ft.get_auxr_data<FPTYPE>()[ir] = std::complex<FPTYPE>(in[ir],0);
        }
        this->ft.fftxyfor(ft.get_auxr_data<FPTYPE>(),ft.get_auxr_data<FPTYPE>());
    }
    this->gatherp_scatters(this->ft.get_auxr_data<FPTYPE>(), this->ft.get_auxg_data<FPTYPE>());
    
    this->ft.fftzfor(ft.get_auxg_data<FPTYPE>(),ft.get_auxg_data<FPTYPE>());

    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += factor / FPTYPE(this->nxyz) * this->ft.get_auxg_data<FPTYPE>()[this->ig2isz[ig]];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = this->ft.get_auxg_data<FPTYPE>()[this->ig2isz[ig]] / FPTYPE(this->nxyz);
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
}

/// 
/// transform reciprocal space to real space
/// f(r)=1/V * \sum_{g} c(k,g)*exp(ig*r)
/// c(k,g)=c_k(g)*exp(ik*r)
/// f(r)=1/V * \sum_{g} c_k(g)*exp(i(g+k)*r)
/// Here we use c(k,g)
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), complex<double>
/// 
template <typename FPTYPE>
void PW_Basis:: recip2real(const std::complex<FPTYPE> * in, std::complex<FPTYPE> * out, const bool add, const FPTYPE factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        ft.get_auxg_data<FPTYPE>()[i] = std::complex<FPTYPE>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.get_auxg_data<FPTYPE>()[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftzbac(ft.get_auxg_data<FPTYPE>(), ft.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->ft.get_auxg_data<FPTYPE>(),this->ft.get_auxr_data<FPTYPE>());

    this->ft.fftxybac(ft.get_auxr_data<FPTYPE>(),ft.get_auxr_data<FPTYPE>());
    
    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] += factor * this->ft.get_auxr_data<FPTYPE>()[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] = this->ft.get_auxr_data<FPTYPE>()[ir];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), double
///
template <typename FPTYPE>
void PW_Basis:: recip2real(const std::complex<FPTYPE> * in, FPTYPE * out, const bool add, const FPTYPE factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        ft.get_auxg_data<FPTYPE>()[i] = std::complex<double>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.get_auxg_data<FPTYPE>()[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftzbac(ft.get_auxg_data<FPTYPE>(), ft.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->ft.get_auxg_data<FPTYPE>(), this->ft.get_auxr_data<FPTYPE>());

    if(this->gamma_only)
    {
        this->ft.fftxyc2r(ft.get_auxr_data<FPTYPE>(),ft.get_rspace_data<FPTYPE>());

        // r2c in place
        const int npy = this->ny * this->nplane;

        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] += factor * this->ft.get_rspace_data<FPTYPE>()[ix*npy + ipy];
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] = this->ft.get_rspace_data<FPTYPE>()[ix*npy + ipy];
                }
            }
        }
    }
    else
    {
        this->ft.fftxybac(ft.get_auxr_data<FPTYPE>(),ft.get_auxr_data<FPTYPE>());
        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] += factor * this->ft.get_auxr_data<FPTYPE>()[ir].real();
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] = this->ft.get_auxr_data<FPTYPE>()[ir].real();
            }
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
}

template void PW_Basis::real2recip<float>(const float * in, std::complex<float> * out, const bool add, const float factor); //in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::real2recip<float>(const std::complex<float> * in, std::complex<float> * out, const bool add, const float factor); //in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::recip2real<float>(const std::complex<float> * in, float *out, const bool add, const float factor); //in:(nz, ns)  ; out(nplane,nx*ny)
template void PW_Basis::recip2real<float>(const std::complex<float> * in, std::complex<float> * out, const bool add, const float factor);

template void PW_Basis::real2recip<double>(const double * in, std::complex<double> * out, const bool add, const double factor); //in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::real2recip<double>(const std::complex<double> * in, std::complex<double> * out, const bool add, const double factor); //in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::recip2real<double>(const std::complex<double> * in, double *out, const bool add, const double factor); //in:(nz, ns)  ; out(nplane,nx*ny)
template void PW_Basis::recip2real<double>(const std::complex<double> * in, std::complex<double> * out, const bool add, const double factor);

}