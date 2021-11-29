#include "pw_basis.h"

#ifdef __MPI
#include "mpi.h"
#include "../src_parallel/parallel_global.h"
#endif
namespace ModulePW
{

//in: (nplane,nxy) out:(nz,nst)
void PW_Basis:: gatherp_scatters(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int bigixy = this->is2bigixy[is];
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                out[is*nz+iz] = in[bigixy*nz+iz];
            }
        }
        return;
    }
#ifdef __MPI
    std::complex<double> * tmp;
    if(this->poolrank == 0) tmp = new std::complex<double> [this->nz * this->nstot];
    
    //gather planes of different processors
    for(int ixy = 0 ; ixy < this->nxy ; ++ixy)
    {
        if(this->ixy2ip[ixy] == -1) continue;
        int istot = 0;
        if(this->poolrank == 0) istot = this->ixy2istot[ixy];
        int bigixy = (ixy / ny) * bigny + ixy % ny;
        MPI_Gatherv(&in[bigixy*this->nplane], this->nplane, mpicomplex, &tmp[istot*this->nz], 
                    this->numz,this->startz,mpicomplex,0,POOL_WORLD);
    }
    
    //scatter sticks to different processors
    MPI_Scatterv(tmp, this->nstnz_per, this->startnsz_per,mpicomplex,out,
                    this->nstnz,mpicomplex,0, POOL_WORLD); 
    
    if(this->poolrank == 0) delete[] tmp;
#endif
    return;
}

//in: (nplane,nxy) out:(nz,nst)
// void PW_Basis:: gatherp_scatters_gamma(complex<double> *in, complex<double> *out)
// {
//     int hx = int ((this->nx + 2)/2);
//     if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
//     {
//         for(int is = 0 ; is < this->nst ; ++is)
//         {
//             int ixy = is2ixy[is];
//             int ix = ixy % this->ny;
//             int iy = int( ixy / this->ny);
//             int ihxy = ix + iy * hx;
//             for(int iz = 0 ; iz < this->nz ; ++iz)
//             {
//                 out[is*nz+iz] = in[ihxy*nz+iz];
//             }
//         }
//         return;
//     }
// #ifdef __MPI
//     std::complex<double> * tmp;
//     if(this->poolrank == 0) tmp = new std::complex<double> [this->nz * this->nstot];
    
//     //gather planes of different processors
//     for(int ixy = 0 ; ixy < this->ny * hx ; ++ixy)
//     {
//         if(this->ixy2ip[ixy] == -1) continue;
//         int istot = 0;
//         if(this->poolrank == 0) istot = this->ixy2istot[ixy];
//         MPI_Gatherv(&in[ixy*this->nplane], this->nplane, mpicomplex, &tmp[istot*this->nz], 
//                     this->numz,this->startz,mpicomplex,0,POOL_WORLD);
//     }
    
//     //scatter sticks to different processors
//     MPI_Scatterv(tmp, this->nstnz_per, this->startnsz_per,mpicomplex,out,
//                     this->nstnz,mpicomplex,0, POOL_WORLD);
//     if(this->poolrank == 0) delete[] tmp;
// #endif
//     return;
// }

//in: (nz,nst) out:(nplane,nxy)
void PW_Basis:: gathers_scatterp(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) //In this case nrxx=nx*ny*nz, nst = nstot, 
    {
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] = 0.0;
        }
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int bigixy = is2bigixy[is];
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                out[bigixy*nz+iz] = in[is*nz+iz];
            }
        }
        return;
    }
#ifdef __MPI
    if(this->poolnproc == 1) return;
    std::complex<double> * tmp;
    if(this->poolrank == 0) tmp = new std::complex<double> [this->nz * this->nstot];

    //scatter sticks to different processors
    MPI_Gatherv(in, this->nstnz, mpicomplex, tmp,
                    this->nstnz_per, this->startnsz_per, mpicomplex, 0, POOL_WORLD);

    for(int ir = 0 ; ir < this->nrxx ; ++ir) out[ir] = 0.0;
    for(int istot = 0 ; istot < this->nstot ; ++istot)
    {
        int bigixy = this->istot2bigixy[istot];
        MPI_Scatterv(&tmp[istot*this->nz], this->numz,this->startz, mpicomplex, &out[bigixy*this->nplane], 
                    this->nplane,mpicomplex,0,POOL_WORLD);
    }
    if(this->poolrank == 0) delete[] tmp;
#endif
    return;
}

//in: (nz,nst) out:(nplane,nxy)
// void PW_Basis:: gathers_scatterp_gamma(complex<double> *in, complex<double> *out)
// {
//     int hx = int ((nx + 2)/2);
//     if(this->poolnproc == 1) //In this case nrxx=nx*ny*nz, nst = nstot, 
//     {
//         for(int ir = 0 ; ir < this->nrxx ; ++ir)
//         {
//             out[ir] = 0.0;
//         }
//         for(int is = 0 ; is < this->nst ; ++is)
//         {
//             int ixy = is2ixy[is];
//             int ix = ixy % this->ny;
//             int iy = int( ixy / this->ny);
//             int ihxy = ix + iy * hx;
//             for(int iz = 0 ; iz < this->nz ; ++iz)
//             {
//                 out[ihxy*nz+iz] = in[is*nz+iz];
//             }
//         }
//         return;
//     }
// #ifdef __MPI
//     if(this->poolnproc == 1) return;
//     std::complex<double> * tmp;
//     if(this->poolrank == 0) tmp = new std::complex<double> [this->nz * this->nstot];

//     //scatter sticks to different processors
//     MPI_Gatherv(in, this->nstnz, mpicomplex, tmp,
//                     this->nstnz_per, this->startnsz_per, mpicomplex, 0, POOL_WORLD);

//     for(int ir = 0 ; ir < this->nrxx ; ++ir) out[ir] = 0.0;
//     for(int istot = 0 ; istot < this->nstot ; ++istot)
//     {
//         int ixy = this->istot2ixy[istot];
//         MPI_Scatterv(&tmp[istot*this->nz], this->numz,this->startz, mpicomplex, &out[ixy*this->nplane], 
//                     this->nplane,mpicomplex,0,POOL_WORLD);
//     }
//     if(this->poolrank == 0) delete []tmp;
// #endif
//     return;
// }

}