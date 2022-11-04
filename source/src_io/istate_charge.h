#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include<vector>
#include<module_base/matrix.h>
#include<module_base/complexmatrix.h>
#include "src_lcao/local_orbital_charge.h"
#include "module_gint/gint_gamma.h"
#include "module_psi/psi.h"

class IState_Charge
{
public:
    IState_Charge(psi::Psi<double>* psi_gamma_in,
        Local_Orbital_Charge &loc_in);
    ~IState_Charge();

	void begin(Gint_Gamma &gg);

private:

	int *bands_picked;

#ifdef __MPI
	void idmatrix(const int &ib);
#endif
    psi::Psi<double>* psi_gamma;
    Local_Orbital_Charge* loc;

};
#endif
