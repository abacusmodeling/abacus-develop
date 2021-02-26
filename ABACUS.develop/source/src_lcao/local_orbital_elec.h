#ifndef LOCAL_ORBTIAL_ELEC
#define LOCAL_ORBITAL_ELEC

#include "../src_pw/tools.h"

class Local_Orbital_Elec
{
	public:

	Local_Orbital_Elec(){};
	~Local_Orbital_Elec(){};

	// mohan add 2021-02-09
	void solve_elec_stru(const int &istep);

	private:

	// set matrix and grid integral
	void set_matrix_grid(void);

	void before_solver(const int &istep);

	void solver(const int &istep); 

};

#endif
