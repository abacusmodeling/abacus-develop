#ifndef DIAGH_H
#define DIAGH_H

#include "module_base/complexmatrix.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"

namespace ModuleHSolver
{

class DiagH
{
    public:

    //virtual void init()=0;

    virtual int diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<std::complex<double>> &phi,
        double *eigenvalue_in)=0;

};

}

#endif