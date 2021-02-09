#include "subgrid_oper.h"
#include "../src_pw/tools.h"
#include "../src_pw/global.h"

SubGrid_oper::SubGrid_oper()
{
	trace_lo_tot = new int[1];	
	lgd=0;
	allocate_totwfc = false;
}
SubGrid_oper::~SubGrid_oper()
{
	delete[] trace_lo_tot;
}

//--------------------------------------------------
// because only the DIAG_WORLD processors
// have augmented wave functions,
// others will not have, we use this 
// augmented wave functions to calculate
// the force, not related to the grid integration
// Only DIAG_WORLD use density matrix (2D) only.
//--------------------------------------------------
