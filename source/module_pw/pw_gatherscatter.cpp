#include "pw_basis.h"

//#ifdef __MPI
#include "mpi.h"
#include "../src_parallel/parallel_global.h"

//in: (nplane,nxy) out:(nz,ns)
void PW_Basis:: gatherp_scatters(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) return;

    complex<double> * tmp;
    if(this->poolrank == 0) tmp = new complex<double> [this->nz * this->nstot];
    
    //gather planes of different processors
    for(int ixy = 0 ; ixy < this->nxy ; ++ixy)
    {
        if(this->ixy2ip[ixy] == -1) continue;
        int istot = 0;
        if(this->poolrank == 0) istot = this->ixy2istot[ixy];
        MPI_Gatherv(&in[ixy*this->nplane], this->nplane, mpicomplex, &tmp[istot*this->nz], 
                    this->numz,this->startz,mpicomplex,0,POOL_WORLD);
    }
    
    //scatter sticks to different processors
    MPI_Scatterv(tmp, this->npw_per, this->startnpw,mpicomplex,out,
                    this->npw,mpicomplex,0, POOL_WORLD); 
    delete[] tmp;
    return;
}

////in: (nz,ns) out:(nplane,nxy)
void gathers_scatterp(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) return;
    complex<double> * tmp;
    if(this->poolrank == 0) tmp = new complex<double> [this->nz * this->nstot];

    //scatter sticks to different processors
    MPI_Gatherv(in, this->npw, mpicomplex, tmp,
                    this->npw_per, this->startnpw, mpicomplex, 0, POOL_WORLD);

    for(int ir = 0 ; ir < this->nrxx ; ++ir) out[ir] = 0.0;
    for(int istot = 0 ; istot < this->nstot ; ++istot)
    {
        int ixy = this->istot2ixy[istot];
        MPI_Scatterv(&tmp[istot*this->nz], this->numz,this->startz, mpicomplex, &out[ixy*this->nplane], 
                    mpicomplex,0,POOL_WORLD);
    }

    return;
}
//#endif