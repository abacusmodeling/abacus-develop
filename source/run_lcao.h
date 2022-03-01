//==========================================================
// AUTHOR : mohan
// DATE : 2021-01-20
//==========================================================
#ifndef RUN_LCAO_H
#define RUN_LCAO_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_orbital/ORB_control.h"
#include "input.h"
#include "module_ensolver/en_solver.h"

class Run_lcao
{

	public:

    Run_lcao();
    ~Run_lcao();

	// perform Linear Combination of Atomic Orbitals (LCAO) calculations
    static void lcao_line(ModuleEnSover::En_Solver *p_ensolver);

private:
    static void Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell_pseudo& ucell);
};

#endif
