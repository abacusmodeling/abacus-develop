#include "sto_wf.h"
#include "global.h"


Stochastic_WF::Stochastic_WF()
{
    chig  = new ComplexMatrix[1];
    chi0  = new ComplexMatrix[1];
}

Stochastic_WF::~Stochastic_WF()
{ 
    delete[] chig;
    delete[] chi0;
}

void Stochastic_WF::init()
{
    //wait for init
    int nrxx;
    int nx,ny,nz;

    //distribute nchi for each process
    nchip = int(nchi/NPROC_IN_POOL);
    if(RANK_IN_POOL < nchi%NPROC_IN_POOL) ++nchip;

    complex<double> ui(0,1);
    
    delete[] chi0;
    //We temporarily init one group of orbitals for all k points.
    //This save memories.
    chi0 = new ComplexMatrix[1]; 
    chi0[0].create(nchip,nrxx,0);
    //init with random number
    for(int i=0; i<chi0[0].size; ++i)
    {
        chi0[0].c[i]=exp(2*PI*rand()*ui);
    }

    delete[] chig;
    int nkk = 1; // We temporarily use gamma k point.
    chig = new ComplexMatrix[1];
    chig[0].create(nchip,nrxx,0);
    return;
}


