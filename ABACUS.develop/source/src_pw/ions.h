#ifndef IONS_H
#define IONS_H

#include "electrons.h"
#include "tools.h"
#include "../src_ions/ions_move_methods.h"
#include "../src_ions/lattice_change_methods.h"
#include "charge_extra.h"
//#include "../src_develop/src_md/md.h"
#include "sto_elec.h" //mohan added 2021-01-28

class Ions
{
public:

    Ions(){};
    ~Ions(){};

    void opt_ions_pw(void);

private:

	// mohan add 2021-01-28
    Electrons elec;

	// mohan add for stochastic wave functions
	Stochastic_Elec elec_sto;	

	// mohan add 2021-01-28
	// mohan moved this variable from electrons.h to ions.h
    int istep;

	Ions_Move_Methods IMM;

	//MD md; //mohan add 2011-11-07

	Charge_Extra CE;
	
	Lattice_Change_Methods LCM;

	bool force_stress(const int &istep, int &force_step, int &stress_step);    // pengfei Li 2018-05-14

	bool force_stress(const int &istep);

    void update_pot(void);

};

#endif
