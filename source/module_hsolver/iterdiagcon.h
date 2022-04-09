#ifndef ITERDIAGCONTROL_H
#define ITERDIAGCONTROL_H

#include "module_psi/psi.h"
#include "module_base/complexmatrix.h"
#include "src_pw/hamilt_pw.h"

namespace ModuleHSolver
{

class IterDiagControl
{
    public: 
    static double PW_DIAG_THR;
    static int PW_DIAG_NMAX;
    /// record for how many bands not have convergence eigenvalues
    static int notconv;
    /// average steps of last cg diagonalization for each band.
    static double avg_iter;
    /// record the times of trying iterative diagonalization
    static int ntry;

    static void diagH_subspace(
        Hamilt_PW* phm,
        const int nstart,
        const int n_band,
        const ModulePsi::Psi<std::complex<double>> &psi,
        ModulePsi::Psi<std::complex<double>> &evc,
        double *en);

    static void diagH_LAPACK(
        const int nstart,
        const int nbands,
        const ModuleBase::ComplexMatrix &hc,
        const ModuleBase::ComplexMatrix &sc,
        const int ldh, // nstart
        double *e,
        ModuleBase::ComplexMatrix &hvec);
};

}

#endif