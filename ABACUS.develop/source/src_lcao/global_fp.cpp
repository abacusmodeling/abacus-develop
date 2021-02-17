#include "global_fp.h"
#include "src_pw/global.h"

Grid_Driver GridD;
Parallel_Orbitals ParaO;
Local_Orbital_Charge LOC;
Local_Orbital_wfc LOWF;
LCAO_Matrix LM;
Use_Hamilt_Matrix UHM;
SubGrid_oper SGO; //mohan add 2012-01-12
Exx_Lcao exx_lcao(exx_global.info); // Peize Lin add 2016-12-03
