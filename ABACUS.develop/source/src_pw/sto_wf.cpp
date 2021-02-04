#include "sto_wf.h"
#include "global.h"


Stochastic_WF::Stochastic_WF()
{
    chi  = new ComplexMatrix[1];
    chi0  = new ComplexMatrix[1];
}

Stochastic_WF::~Stochastic_WF()
{ 
    delete[] chi;
    delete[] chi0;
}

void Stochastic_WF::init()
{
    int nrxx = 1;
    complex<double> ui(0,1);
    
    delete[] chi0;
    //We temporarily init one group of orbitals for all k points.
    //This save memories.
    chi0 = new ComplexMatrix[1]; 
    chi0[0].create(nchi,nrxx,0);
    //init with random number
    for(int i=0; i<chi0[0].size; i++)
    {
        chi0[0].c[i]=exp(2*PI*rand()*ui);
    }

    delete[] chi;
    int nkk = 1; // We temporarily use gamma k point.
    chi = new ComplexMatrix[nkk];
    return;
}



