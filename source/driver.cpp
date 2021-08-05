#include "driver.h"
#include "run_pw.h"
#include "input.h"
#include "input_conv.h"
#ifdef __LCAO
#include "run_lcao.h"
#include "src_lcao/global_fp.h"
#endif
#include "src_pw/global.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"

Driver::Driver(){}


Driver::~Driver(){}


void Driver::init()
{
	TITLE("Driver","init");

	// (1) read the input parameters.
	this->reading();

	// (2) welcome to the std::atomic world!
	this->atomic_world();

	// (3) close all of the running logs 
	INPUT.close_log();

	return;
}


void Driver::reading(void)
{
	timer::tick("Driver","reading");

	// (1) read INPUT 
	INPUT.Init( GlobalV::global_in_card );

	// (2) copy the variables from INPUT to each class
	Input_Conv::Convert();

	// (3) define the 'DIAGONALIZATION' world in MPI
	Parallel_Global::split_diag_world(GlobalV::DIAGO_PROC);
	Parallel_Global::split_grid_world(GlobalV::DIAGO_PROC);
	OUT(GlobalV::ofs_running,"DRANK",GlobalV::DRANK+1);
	OUT(GlobalV::ofs_running,"DSIZE",GlobalV::DSIZE);
	OUT(GlobalV::ofs_running,"DCOLOR",GlobalV::DCOLOR+1);
	OUT(GlobalV::ofs_running,"GRANK",GlobalV::GRANK+1);
	OUT(GlobalV::ofs_running,"GSIZE",GlobalV::GSIZE);

#ifdef __MPI
    // (4)  divide the GlobalV::NPROC processors into GlobalV::NPOOL for k-points parallelization.
    GlobalC::Pkpoints.init_pools();
#endif

    // (5) Read in parameters about wannier functions.
    winput::Init( GlobalV::global_wannier_card );

    // (6) Print the parameters into INPUT file.
    std::stringstream ss1;
    ss1 << GlobalV::global_out_dir << GlobalV::global_in_card;
    INPUT.Print( ss1.str() );
    //DONE(GlobalV::ofs_running,"READING CARDS");

	timer::tick("Driver","reading");
	return;
}

void Driver::atomic_world(void)
{
	TITLE("Driver","atomic_world");

	//--------------------------------------------------
	// choose basis sets:
	// pw: plane wave basis std::set
	// lcao_in_pw: LCAO expaned by plane wave basis std::set
	// lcao: linear combination of std::atomic orbitals
	//--------------------------------------------------
	if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		Run_pw::plane_wave_line();
	}
#ifdef __LCAO
	else if(GlobalV::BASIS_TYPE=="lcao")
	{
		Run_lcao::lcao_line();
	}
#endif

	timer::finish( GlobalV::ofs_running );

	Memory::print_all( GlobalV::ofs_running ) ;

	return;
}
