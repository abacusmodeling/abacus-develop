#ifndef RUN_MD_LCAO_H
#define RUN_MD_LCAO_H 

#include "../src_lcao/LOOP_elec.h"
#include "../src_ions/ions_move_methods.h"
#include "../src_pw/charge_extra.h"
#include "MD_basic.h"
#include "../src_ions/lattice_change_methods.h"

class Run_MD_LCAO
{

	public:

	Run_MD_LCAO();
	~Run_MD_LCAO();

	LOOP_elec LOE;

	void opt_cell(void);
	void opt_ions(void);

	private:

	Ions_Move_Methods IMM;

	//bool force_stress(void);
	Lattice_Change_Methods LCM;

	int istep;

	// electron charge density extropolation method
	Charge_Extra CE;

	void final_scf(void);

};

#endif
