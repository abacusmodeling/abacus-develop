#ifndef IONS_H
#define IONS_H

#include "electrons.h"
#include "tools.h"
#include "../src_ions/ions_move_methods.h"
#include "../src_ions/lattice_change_methods.h"
#include "charge_extra.h"
//#include "../src_develop/src_md/md.h"

class Ions : public electrons
{
public:

    Ions(){};
    ~Ions(){};

    void opt_ions_pw(void);

private:

	Ions_Move_Methods IMM;
	//MD md; //mohan add 2011-11-07
	Charge_Extra CE;
	
	Lattice_Change_Methods LCM;
	bool force_stress(const int &istep, int &force_step, int &stress_step);    // pengfei Li 2018-05-14

	bool force_stress(const int &istep);
    void update_pot(void);
    void extrapolate_wfcs(void);

};

#endif
