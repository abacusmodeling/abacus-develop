#ifndef DIAGH_H
#define DIAGH_H

#include "module_base/complexmatrix.h"

namespace ModuleHSolver
{

class DiagH
{
    public:

    virtual void init()=0;

    virtual int diag(
        const int &dim_in,
        const double *precondition_in,
        ModuleBase::ComplexMatrix &phi,
        double *eigenvalue_in)=0;

    void schmit_orth
    (
        const int& dim,
        const int& m,     //end
        const ModuleBase::ComplexMatrix &psi,
        std::complex<double> *sphi,
        std::complex<double> *psi_m
    );
};

}

#endif