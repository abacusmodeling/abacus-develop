#include "pw_basis.h"
//
//distribute real-space grids to different processors
//Known: nx, ny, nz, poolnproc, poolrank
//output: nrxx, startz, numz
//
void PW_Basis::distribute_r()
{
    int npz = this->nz / this->poolnproc;
    int modz = this->nz % this->poolnproc;
    this->startz[0] = 0;
    for(int ip = 0 ; ip < this->poolnproc ; ++ip)
    {
        this->numz[ip] = npz;
        if(ip < modz)   this->numz[ip]++;
        if(ip < this->poolnproc - 1)   this->startz[ip+1] += numz[ip];
    }
    this->nrxx = this->numz[this->poolrank];
    return;
}

void PW_Basis::gatherplane()
{
    
}