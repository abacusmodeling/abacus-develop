#include "./stress_func.h"
#include "./xc_functional.h"
#include "./myfunc.h"
#include "./xc_gga_pw.h"

//calculate the mGGA stress correction in PW and LCAO
void Stress_Func::stress_mgga(matrix& sigma) 
{
	timer::tick("Stress_Func","stress_mgga");

	int current_spin = 1;
	if (GlobalV::NSPIN==4) WARNING_QUIT("stress_mgga","noncollinear stress + mGGA not implemente");

	timer::tick("Stress_Func","stress_mgga");
	return;
}
