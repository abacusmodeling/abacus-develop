#ifndef LOOP_CELL_H
#define LOOP_CELL_H
#include<vector>
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_orbital/ORB_control.h"
#include "src_lcao/LCAO_matrix.h"
#include "module_ensolver/en_solver.h"

class LOOP_cell
{
	public:

	LOOP_cell(Parallel_Orbitals &pv);
	~LOOP_cell();

	void opt_cell(ORB_control &orb_con, ModuleEnSover::En_Solver *p_ensolver);

private:
    LCAO_Matrix LM;
};

#endif
