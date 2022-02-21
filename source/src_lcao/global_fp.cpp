#include "global_fp.h"
#include "../src_pw/global.h"

namespace GlobalC
{
Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver,GlobalV::test_grid);
Pdiag_Double ParaO;
LCAO_Matrix LM;
SubGrid_oper SGO; //mohan add 2012-01-12
Exx_Lcao exx_lcao(GlobalC::exx_global.info); // Peize Lin add 2016-12-03
}
