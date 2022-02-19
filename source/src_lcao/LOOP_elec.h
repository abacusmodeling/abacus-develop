#ifndef LOOP_ELEC_H
#define LOOP_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"

class LOOP_elec
{
	public:

	LOOP_elec(){};
	~LOOP_elec(){};

	// mohan add 2021-02-09
    void solve_elec_stru(const int& istep,
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::matrix>& dm_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        std::vector<ModuleBase::ComplexMatrix>& dm_k);

	private:

	// set matrix and grid integral
	void set_matrix_grid(void);

    void before_solver(const int& istep,
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k);

    void solver(const int& istep,
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::matrix>& dm_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        std::vector<ModuleBase::ComplexMatrix>& dm_k);

};

#endif
