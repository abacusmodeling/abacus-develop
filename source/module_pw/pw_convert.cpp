#include "fft.h"
#include <complex>
#include "pw_basis.h"
#include "../module_base/global_function.h"

//
//transform real space to reciprocal space
//in: (nplane,nx,ny)
//out: (nz, ns)
//
void PW_Basis:: real2recip(complex<double> * in, complex<double> * out)
{
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.c_rspace[ir] = in[ir];
    }
    this->ft.executefftw("2for");

    this->gatherp_scatters(this->ft.c_rspace, this->ft.c_gspace);
    
    this->ft.executefftw("1for");

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.c_gspace[this->ig2fft[ig]];
    }
    return;
}

//
//transform real space to reciprocal space
//in: (nplane,nx,ny)
//out: (nz, ns)
//
void PW_Basis:: real2recip(double * in, complex<double> * out)
{
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.r_rspace[ir] = in[ir];
    }
    this->ft.executefftw("2r2c");

    this->gatherp_scatters(this->ft.c_rspace, this->ft.c_gspace)
    
    this->ft.executefftw("1for");

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.c_gspace[this->ig2fft[ig]];
    }
    return;
}

//
//transform real space to reciprocal space
//in: (nz,ns)
//out: (nplane, nx, ny)
//
void PW_Basis:: recip2real(complex<double> * in, complex<double> * out)
{
    for(int igg = 0 ; igg < this->ns * this->nz ; ++igg)
    {
        this->ft.c_gspace[igg] = 0.0;
    }
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.c_gspace[this->ig2fft[ig]] = in[ig];
    }
    this->ft.executefftw("1bac");

    this->gathers_scatterp(this->ft.c_gspace,this->ft.c_rspace);
    
    this->ft.executefftw("2bac");

    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.c_rspace[ir];
    }
    return;
}

//
//transform real space to reciprocal space
//in: (nz,ns)
//out: (nplane, nx, ny)
//
void PW_Basis:: recip2real(complex<double> * in, double * out)
{
    for(int igg = 0 ; igg < this->ns * this->nz ; ++igg)
    {
        this->ft.c_gspace[igg] = 0.0;
    }
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.c_gspace[this->ig2fft[ig]] = in[ig];
    }
    this->ft.executefftw("1bac");
    
    this->gathers_scatterp(this->ft.c_gspace, this->ft.c_rspace)

    this->ft.executefftw("2c2r");

    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.r_rspace[ir];
    }
    return;
}