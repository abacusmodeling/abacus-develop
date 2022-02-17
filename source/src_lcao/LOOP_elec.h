#ifndef LOOP_ELEC_H
#define LOOP_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "wfc_dm_2d.h"

class LOOP_elec
{
	public:

	LOOP_elec(){};
	~LOOP_elec(){};

	// mohan add 2021-02-09
    void solve_elec_stru(const int& istep,
        Wfc_Dm_2d& wfc_dm_2d);

	private:

	// set matrix and grid integral
	void set_matrix_grid(void);

	void before_solver(const int &istep, Wfc_Dm_2d &wfc_dm_2d);

    void solver(const int& istep,
        Wfc_Dm_2d &wfc_dm_2d);

};

#endif
