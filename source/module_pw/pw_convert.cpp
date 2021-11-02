#include "fft.h"
#include <complex>
#include "pw_basis.h"
#include "../module_base/global_function.h"

void PW_Basis:: real2recip(complex<double> * in, complex<double> * out)
{
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.c_rspace[ir] = in[ir];
    }
    this->ft.executefftw("2for");

    for(int is = 0 ; is < this->ns ; ++is)
    {
        int ixy = is2ixy[is];
        for(int iz = 0 ; iz < this->nz ; ++iz)
        {
            this->ft.c_gspace[is*nz+iz] = this->ft.c_rspace[ixy*nz+iz];
        }
    }
    this->ft.executefftw("1for");

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.c_gspace[this->ig2fft[ig]];
    }
}

void PW_Basis:: real2recip(double * in, complex<double> * out)
{
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.r_rspace[ir] = in[ir];
    }
    this->ft.executefftw("2r2c");

    int hx = int ((nx + 2)/2);
    for(int is = 0 ; is < this->ns ; ++is)
    {
        int ixy = is2ixy[is];
        int ix = ixy % this->ny;
        int iy = int( ixy / this->ny);
        int ihxy = ix + iy * hx;
        for(int iz = 0 ; iz < this->nz ; ++iz)
        {
            this->ft.c_gspace[is*nz+iz] = this->ft.c_rspace[ihxy*nz+iz];
        }
    }
    this->ft.executefftw("1for");

    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        out[ig] = this->ft.c_gspace[this->ig2fft[ig]];
    }
}

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
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.c_rspace[ir] = 0.0;
    }
    for(int is = 0 ; is < this->ns ; ++is)
    {
        int ixy = is2ixy[is];
        for(int iz = 0 ; iz < this->nz ; ++iz)
        {
            this->ft.c_rspace[ixy*nz+iz] = this->ft.c_gspace[is*nz+iz];
        }
    }
    this->ft.executefftw("2bac");

    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.c_rspace[ir];
    }
}

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
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.c_rspace[ir] = 0.0;
    }
    int hx = int ((nx + 2)/2);
    for(int is = 0 ; is < this->ns ; ++is)
    {
        int ixy = is2ixy[is];
        int ix = ixy % this->ny;
        int iy = int( ixy / this->ny);
        int ihxy = ix + iy * hx;
        for(int iz = 0 ; iz < this->nz ; ++iz)
        {
            this->ft.c_rspace[ihxy*nz+iz] = this->ft.c_gspace[is*nz+iz];
        }
    }
    this->ft.executefftw("2c2r");

    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.r_rspace[ir];
    }
}