#include "LCAO_deepks_test.h"
#ifdef __MPI
#include <mpi.h>
#endif

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Test of module_deepks" << std::endl;

#ifdef __MPI
	MPI_Init(&argc,&argv);
#endif
    calculate();
#ifdef __MPI
	MPI_Finalize();
#endif

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

	std::cout << "Total checks : " << test.total_check <<std::endl;
	std::cout << "Failed checks : " << test.failed_check <<std::endl;

    return;
}