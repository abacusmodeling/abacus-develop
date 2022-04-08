#ifndef ITERDIAGCONTROL_H
#define ITERDIAGCONTROL_H

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

};

}

#endif