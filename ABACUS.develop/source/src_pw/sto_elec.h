#ifndef INCLUDE_STO_ELEC_H
#define INCLUDE_STO_ELEC_H

#include "tools.h"
#include "threshold_elec.h"
#include "sto_wf.h"

//----------------------------------------------
// methods based on stochastic wave functions
//----------------------------------------------

class Stochastic_Elec: private Threshold_Elec
{

	public:

    // constructor and deconstructor
    Stochastic_Elec();
    ~Stochastic_Elec();

    int iter;

    static double avg_iter;

    void scf_stochastic(const int &istep);

	private:

	Stochastic_WF swf;

    void c_bands(const int &istep);

    bool check_stop_now(void);
};

#endif// Eelectrons_Stochastic
