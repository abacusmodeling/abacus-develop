#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include<vector>
#include<module_base/matrix.h>
#include<module_base/complexmatrix.h>

class IState_Charge
{
public:
    IState_Charge(std::vector<ModuleBase::matrix> *wfc_gamma_in,
        std::vector<ModuleBase::matrix> *dm_gamma_in);
    ~IState_Charge();

	void begin();

private:

	int *bands_picked;

#ifdef __MPI
	void idmatrix(const int &ib);
#endif
    std::vector<ModuleBase::matrix> *wfc_gamma;
    std::vector<ModuleBase::matrix> *dm_gamma;

};
#endif
