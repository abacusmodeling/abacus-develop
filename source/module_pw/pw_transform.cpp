#include "fft.h"
#include <complex>
#include "pw_basis.h"
#include <cassert>
#include "../module_base/global_function.h"

namespace ModulePW
{

//
//transform real space to reciprocal space
//in: (nplane,nx,ny)
//out: (nz, ns)
//
void PW_Basis:: real2recip(std::complex<double> * in, std::complex<double> * out)
{
    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.c_gspace[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.c_gspace,ft.c_gspace);

    this->gatherp_scatters(this->ft.c_gspace, this->ft.c_gspace2);
    
    this->ft.fftzfor(ft.c_gspace2,ft.c_gspace);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.c_gspace[this->ig2isz[ig]];
    }
    return;
}

//
//transform real space to reciprocal space
//in: (nplane,nx,ny)
//out: (nz, ns)
//
void PW_Basis:: real2recip(double * in, std::complex<double> * out)
{
    assert(this->gamma_only == true);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.r_rspace[ir] = in[ir];
    }
    // r2c in place
    // int padny = 2 * ny;
    // int nyp = this->bigny * this->nplane;
    // for(int ix = 0 ; ix < this->nx ; ++ix)
    // {
    //     for(int iy = 0 ; iy < this->bigny ; ++iy)
    //     {
    //         for(int iz = 0 ; iz < this->nplane; ++iz)
    //         {
    //             this->ft.r_rspace[ix*padny * nplane * 2 + iy * this->nplane * 2 + iz * 2 ] = in[ix*nyp + iy * this->nplane + iz];
    //         }
    //     }
    // }

    this->ft.fftxyr2c(ft.r_rspace,ft.c_gspace);

    this->gatherp_scatters(this->ft.c_gspace, this->ft.c_gspace2);
    
    this->ft.fftzfor(ft.c_gspace2,ft.c_gspace);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.c_gspace[this->ig2isz[ig]];
    }
    return;
}

//
//transform real space to reciprocal space
//in: (nz,ns)
//out: (nplane, nx, ny)
//
void PW_Basis:: recip2real(std::complex<double> * in, std::complex<double> * out)
{
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.c_gspace, this->nst * this->nz);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.c_gspace[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftzbac(ft.c_gspace, ft.c_gspace2);

    this->gathers_scatterp(this->ft.c_gspace2,this->ft.c_gspace);

    this->ft.fftxybac(ft.c_gspace,ft.c_gspace);
    
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.c_gspace[ir] / double(this->bignxyz);
    }

    return;
}

//
//transform real space to reciprocal space
//in: (nz,ns)
//out: (nplane, nx, ny)
//
void PW_Basis:: recip2real(std::complex<double> * in, double * out)
{
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.c_gspace, this->nst * this->nz);

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.c_gspace[this->ig2isz[ig]] = in[ig];
    }
   this->ft.fftzbac(ft.c_gspace, ft.c_gspace2);
    
    this->gathers_scatterp(this->ft.c_gspace2, this->ft.c_gspace);

    this->ft.fftxyc2r(ft.c_gspace,ft.r_rspace);

    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.r_rspace[ir] / this->bignxyz;
    }
    // c2r in place
    // int padny = 2 * ny;
    // int nyp = this->bigny * this->nplane;
    // for(int ix = 0 ; ix < this->nx ; ++ix)
    // {
    //     for(int iy = 0 ; iy < this->bigny ; ++iy)
    //     {
    //         for(int iz = 0 ; iz < this->nplane; ++iz)
    //         {
    //             out[ix*nyp + iy * this->nplane + iz] = this->ft.r_rspace[ix*padny * nplane *2 + iy * this->nplane*2 + iz*2] / double(this->bignxyz);
    //         }
    //     }
    // }
    return;
}

}