#ifndef INCLUDE_ELECTRONS_STOCHASTIC_H
#define INCLUDE_ELECTRONS_STOCHASTIC_H

#include "tools.h"
#include "threshold_elec.h"

//----------------------------------------------
// methods based on stochastic wave functions
//----------------------------------------------

class Electrons_Stochastic: private Threshold_Elec
{

	public:

    // constructor and deconstructor
    Electrons_Stochastic();
    ~Electrons_Stochastic();

    int iter;
    static double avg_iter;
    int test;

    void scf_stochastic(const int &istep);

	private:

    void c_bands(const int &istep);

    bool check_stop_now(void);
};

#endif// Eelectrons_Stochastic
