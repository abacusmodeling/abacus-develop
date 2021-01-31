#ifndef INCLUDE_ELECTRONS_H
#define INCLUDE_ELECTRONS_H

#include "tools.h"
#include "threshold_elec.h"

class electrons: private Threshold_Elec
{

public:

    // constructor and deconstructor
    electrons();
    ~electrons();

public:

    //==========================================================
    // lberry   : If TRUE Berry phase polarization is calculated.
    // iter     : Current step of electrons iteration.
    // ntry     :
    // notconv  :
    // ethr		: the convergence threshold for eigenvalues
    // avg_iter : Average iteration number in ccgdiagg
    //==========================================================

    int iter;
    static double avg_iter;
    int test;
    int unit;
    double delta_total_energy;

protected:
    void self_consistent(const int &istep);
    void non_self_consistent(void);
    int istep;

private:

    void c_bands(void);
    void compute_magnetization(void);

private:

    // used in electrons()
    bool check_stop_now(void);
    void init_mixstep_final_scf(void);
};

#endif
