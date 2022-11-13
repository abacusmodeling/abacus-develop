#include "global_fp.h"
#include "../src_pw/global.h"

namespace GlobalC
{
Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver,GlobalV::test_grid);
SubGrid_oper SGO; //mohan add 2012-01-12

#ifdef __EXX
Exx_Lcao exx_lcao(GlobalC::exx_info.info_global); // Peize Lin add 2016-12-03
Exx_LRI<double> exx_lri_double(GlobalC::exx_info.info_ri); // Peize Lin add 2022-08-06
Exx_LRI<std::complex<double>> exx_lri_complex(GlobalC::exx_info.info_ri); // Peize Lin add 2022-08-06
#endif
}
