#include "global_fp.h"
#include "../src_pw/global.h"

namespace GlobalC
{
Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver,GlobalV::test_grid);
Parallel_Orbitals ParaO;
Local_Orbital_Charge LOC;
Local_Orbital_wfc LOWF;
LCAO_Matrix LM;
LCAO_Hamilt UHM;
SubGrid_oper SGO; //mohan add 2012-01-12
}
Exx_Lcao exx_lcao(GlobalC::exx_global.info); // Peize Lin add 2016-12-03
