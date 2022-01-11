#include "LCAO_deepks_test.h"

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Hello, this is the DEEPKS module of ABACUS." << std::endl;

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

	std::cout << "--------------------" << std::endl;
	std::cout << " Have a great day! " << std::endl;
	std::cout << "--------------------" << std::endl;

    return;
}