#ifndef LOOP_ELEC_H
#define LOOP_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"

class LOOP_elec
{
	public:

	LOOP_elec(){};
	~LOOP_elec(){};

	// mohan add 2021-02-09
	void solve_elec_stru(const int &istep);

	private:

	// set matrix and grid integral
	void set_matrix_grid(void);

	void before_solver(const int &istep);

	void solver(const int &istep); 

};

#endif
