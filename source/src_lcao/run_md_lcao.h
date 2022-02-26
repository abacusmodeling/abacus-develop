#ifndef RUN_MD_LCAO_H
#define RUN_MD_LCAO_H 

#include "../src_pw/charge_extra.h"
#include "module_orbital/ORB_control.h"
#include "src_lcao/LCAO_matrix.h"
#include "module_ensolver/en_solver.h"

class Run_MD_LCAO
{

	public:

	Run_MD_LCAO(Parallel_Orbitals &pv);
	~Run_MD_LCAO();

	void opt_cell(ORB_control &orb_con, ModuleEnSover::En_Solver *p_ensolver);
	void opt_ions(ModuleEnSover::En_Solver *p_ensolver);
	void md_force_virial(ModuleEnSover::En_Solver *p_ensolver,
		const int &istep,
        const int& numIon, 
        double &potential, 
        ModuleBase::Vector3<double>* force, 
        ModuleBase::matrix& virial);

	private:

	// electron charge density extropolation method
	Charge_Extra CE;
    bool cellchange;
    LCAO_Matrix LM_md;
};

#endif
