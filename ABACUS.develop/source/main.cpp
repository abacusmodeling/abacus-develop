//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================
#include "dc_driv.h"
#include "src_parallel/parallel_global.h"
#include "src_global/global_variable.h"
#include <ctime>
void calculate();

int main(int argc, char **argv)
{
    Parallel_Global::read_mpi_parameters(argc,argv);
    calculate();
    return 0;
}

void calculate()
{
	// divide and conqure (DC) starts here,
	DC_Driv DD;
	DD.init();

#ifdef __MPI
    MPI_Finalize();
#endif
    return;
}
