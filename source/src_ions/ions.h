#ifndef IONS_H
#define IONS_H

#include "../src_pw/electrons.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../src_pw/charge_extra.h"
#include "ions_move_methods.h"
#include "lattice_change_methods.h"
#include "module_esolver/esolver.h"

class Ions
{

	public:

    Ions(){};
    ~Ions(){};

    void opt_ions(ModuleESolver::ESolver *p_esolver);

	private:

	// mohan add 2021-01-28
    Electrons elec;

	// mohan add 2021-01-28
	// mohan moved this variable from electrons.h to ions.h
    int istep;

	Ions_Move_Methods IMM;

	//MD md; //mohan add 2011-11-07

	Charge_Extra CE;
	
	Lattice_Change_Methods LCM;

	//seperate force_stress function first
	bool after_scf(ModuleESolver::ESolver *p_esolver,const int &istep, int &force_step, int &stress_step);
	bool if_do_relax();
	bool if_do_cellrelax();
	bool do_relax(const int& istep, int& jstep, const ModuleBase::matrix& ionic_force, const double& total_energy);
	bool do_cellrelax(const int& istep, const ModuleBase::matrix& stress, const double& total_energy);
	void reset_after_relax(const int& istep);
	void reset_after_cellrelax(int& force_step, int& stress_step, ModuleESolver::ESolver *p_esolver);

    void update_pot(void);

};

#endif
