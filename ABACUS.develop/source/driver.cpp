#include "driver.h"
#include "run_pw.h"
#include "run_lcao.h"
#include "input.h"
#include "input_conv.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"
#include "src_io/print_info.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"

Driver::Driver(){}


Driver::~Driver(){}


void Driver::init()
{
	TITLE("Driver","init");

	// (1) read the input parameters.
	this->reading();


	// (2) welcome to the atomic world!
	this->atomic_world();


	// (3) close all of the running logs 
	INPUT.close_log();

	return;
}


void Driver::reading(void)
{
	timer::tick("Driver","reading",'A');

//---------------------------------------------------------------------------------

	// (1) read INPUT 
	INPUT.Init( global_in_card );

	// (2) copy the variables from INPUT to each class
	Input_Conv::Convert();

	// (3) define the 'DIAGONALIZATION' world in MPI
	Parallel_Global::split_diag_world(DIAGO_PROC);
	Parallel_Global::split_grid_world(DIAGO_PROC);
	OUT(ofs_running,"DRANK",DRANK+1);
	OUT(ofs_running,"DSIZE",DSIZE);
	OUT(ofs_running,"DCOLOR",DCOLOR+1);
	OUT(ofs_running,"GRANK",GRANK+1);
	OUT(ofs_running,"GSIZE",GSIZE);

#ifdef __MPI
    // (4)  divide the NPROC processors into NPOOL for k-points parallelization.
    Pkpoints.init_pools();
#endif

    // (5) Read in parameters about wannier functions.
    winput::Init( global_wannier_card );

    // (6) Print the parameters into INPUT file.
    stringstream ss1;
    ss1 << global_out_dir << global_in_card;
    INPUT.Print( ss1.str() );
    //DONE(ofs_running,"READING CARDS");

    // (7) Setup the unitcell.
	// improvement: a) separating the first reading of the atom_card and subsequent
	// cell relaxation. b) put NLOCAL and NBANDS as input parameters
    ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running);
    //ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running, NLOCAL, NBANDS);
    DONE(ofs_running, "SETUP UNITCELL");

    // (8) Symmetry analysis. 
	// symmetry analysis should be performed every time the cell is changed
    if (SYMMETRY)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }
    
	// (9) Setup the k points according to symmetry.
	kv.set( symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");

	// (10) print information
	// mohan add 2021-01-30
	Print_Info PI;
	PI.setup_parameters();

	timer::tick("Driver","reading",'A');
	return;
}




void Driver::atomic_world(void)
{
	TITLE("Driver","atomic_world");
	timer::tick("Driver","atomic_world",'A');

	//--------------------------------------------------
	// choose basis sets:
	// pw: plane wave basis set
	// lcao_in_pw: LCAO expaned by plane wave basis set
	// lcao: linear combination of atomic orbitals
	//--------------------------------------------------
	if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw")
	{
		Run_pw::plane_wave_line();
	}
	else if(BASIS_TYPE=="lcao")
	{
		Run_lcao::lcao_line();
	}

	timer::tick("Driver","atomic_world",'A');
	timer::finish( ofs_running );

	Memory::print_all( ofs_running ) ;

	return;
}



