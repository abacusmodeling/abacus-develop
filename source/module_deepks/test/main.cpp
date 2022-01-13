#include "LCAO_deepks_test.h"

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Test of module_deepks" << std::endl;

    calculate();

    return 0;
}

void calculate()
{
	test_deepks test;

	test.preparation();

	test.check_dstable();
	test.check_psialpha();

	test.check_pdm();
	test.check_gdmx();

	test.check_descriptor();
	test.check_gvx();

	test.check_edelta();
	test.check_v_delta();
	test.check_e_deltabands();
	test.check_f_delta();

    return;
}