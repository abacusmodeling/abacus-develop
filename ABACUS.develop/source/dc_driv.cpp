#include "dc_driv.h"
#include "run_pw.h"
#include "run_lcao.h"
#include "input.h"
#include "input_conv.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"
#include "src_global/print_info.h"

#include "src_pw/cal_test.h"
#include "src_pw/cal_test0.h"
#include "src_pw/winput.h"

DC_Driv::DC_Driv(){}


DC_Driv::~DC_Driv(){}


void DC_Driv::init()
{
	TITLE("DC_Driv","init");

	time_t time_start = std::time(NULL);

	timer::start();


	// (1) read the input parameters.
	this->reading();


	// (2) welcome to the atomic world!
	this->atomic_world();


	time_t	time_finish= std::time(NULL);

	// print out information before ABACUS ends
	cout << "\n START  Time  : " << ctime(&time_start);
	cout << " FINISH Time  : " << ctime(&time_finish);
	cout << " TOTAL  Time  : " << difftime(time_finish, time_start) << endl;
	cout << " SEE INFORMATION IN : "<<global_out_dir<<endl;

	ofs_running << "\n Start  Time  : " << ctime(&time_start);
	ofs_running << " Finish Time  : " << ctime(&time_finish);

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	ofs_running << " Total  Time  : " << hour << " h "
	            << mins << " mins "
	            << secs << " secs "<< endl;

	INPUT.close_log();

	return;
}


void DC_Driv::reading(void)
{
	timer::tick("DC_Driv","reading",'A');

	// reading input files from the 'INPUT' object
	// the 'INPUT'  is global and can be used anywhere,
	// although I suggest you keep the parameters to be as
	// local as possible -- mohan 2021-01-31
	INPUT.Init( global_in_card );

	// copy the variables from INPUT to each class
	Input_Conv::Convert();
	Input_Conv::Convert_FP();

	// there is a 'DIAGONALIZATION' world I define when I was young
	// the 'DIAGO' world only uses computational resources that 
	// are needed for diagonalization of the Hamiltionian -- mohan 2021-01-31
	Parallel_Global::split_diag_world(DIAGO_PROC);
	Parallel_Global::split_grid_world(DIAGO_PROC);
	OUT(ofs_running,"DRANK",DRANK+1);
	OUT(ofs_running,"DSIZE",DSIZE);
	OUT(ofs_running,"DCOLOR",DCOLOR+1);
	OUT(ofs_running,"GRANK",GRANK+1);
	OUT(ofs_running,"GSIZE",GSIZE);

	//----------------------------------
	// call frag_init
	// * divide the k-points into NPOOL
	// * setup the unit cell
	// * do symmetry analysis
	// * setup the k-points
	//----------------------------------
	this->prepare();


	// for LCAO basis, reading the orbitals and construct
	// the interpolation tables.
	// this part should be moved somewher else -- mohan 2021-01-30
	if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-01
	{
		// read orbital information.
		// init overlap matrix table, which is 'S Table'
		// init kinetical matrix element table, which is 'T Table'
		// init non-local pseudopotential matrix element table, which is 'NL Table'
		hm.hon.set_orb_tables();

		// xiaohui add 2015-09-06
		// (1) divide the H and S matrix into each CPU, count the dimensions
		// (2) set the 'trace' between local H/S and global H/S 
		// (2) allocate the needed H and S memory
		LM.divide_HS_in_frag(); 
	}


	// mohan add 2021-01-30
	Print_Info PI;
	PI.screen_output();



	// Initalize the plane wave basis set
	pw.gen_pw(ofs_running, ucell, kv);
	DONE(ofs_running,"INIT PLANEWAVE");
	cout << " UNIFORM GRID DIM     : " << pw.nx <<" * " << pw.ny <<" * "<< pw.nz << endl;
	cout << " UNIFORM GRID DIM(BIG): " << pw.nbx <<" * " << pw.nby <<" * "<< pw.nbz << endl;


	// mohan add 2010-10-10, just to test the symmetry of a variety
	// of systems.
	//xiaohui modified 2013-03-23,adding "/*"
	if(CALCULATION == "test")
	{
		Cal_Test::adjacent_atoms();
		//Cal_Test::sparsity();
		Cal_Test::test_memory();
		QUIT();
	}

	// mohan add 2010-09-13
	// init the grid, then the charge
	// on grid can be distributed.
	Pgrid.init(pw.ncx, pw.ncy, pw.ncz, pw.nczp, pw.nrxx, pw.nbz, pw.bz); // mohan add 2010-07-22, update 2011-05-04
	

	timer::tick("DC_Driv","reading",'A');
	return;
}




void DC_Driv::atomic_world(void)
{
	TITLE("DC_Driv","atomic_world");
	timer::tick("DC_Driv","atomic_world",'A');

	//xiaohui add 2013-09-01
	if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw")
	{
		Run_pw rp;
		rp.plane_wave_line();
	}
	else if(BASIS_TYPE=="lcao")
	{
		Run_lcao rl;
		rl.lcao_line();
	}

	timer::tick("DC_Driv","atomic_world",'A');
	timer::finish( ofs_running );

	Memory::print_all( ofs_running ) ;

	return;
}


// * initialize the NPOOL for k-points
// * initialize the input parameters for wannier functions
// * setup the unit cell
// * symmetry analysis
// * setup the k-points according to symmetry annalysis
void DC_Driv::prepare(void)
{
	TITLE("Run_Frag","init");
	timer::tick("Run_Frag","init",'B');

#ifdef __MPI
    // (1) If k point number > 1, After reading in NPOOL, 
	// divide the NPROC processprs into NPOOL.
    Pkpoints.init();
#endif

    // (2) Read in parameters about wannier functions.
    winput::Init( global_wannier_card );

    //xiaohui move 3 lines, 2015-09-30
    //stringstream ss2;
    //ss2 << global_out_dir << "INPUTw";
    //winput::Print( ss2.str() );

    // (3) Print the parameters into INPUT file.
    stringstream ss1;
    ss1 << global_out_dir << global_in_card;
    INPUT.Print( ss1.str() );
    //DONE(ofs_running,"READING CARDS");

    // (4) Setup the unitcell.
    ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running);
    DONE(ofs_running, "SETUP UNITCELL");

    // (5) symmetry analysize.
    if (SYMMETRY)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }
    
	// (6) Setup the k points according to symmetry.
	kv.set( symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");
	

	// (7) check the number of basis, the warning should be moved to 
	// other places -- mohan 2021-01-30
	// mohan add 2011-01-5
	if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") 
	{
		if( NLOCAL < NBANDS )
		{
			WARNING_QUIT("UnitCell_pseudo::cal_nwfc","NLOCAL < NBANDS");
		}
		else
		{
			//OUT(ofs_running,"NLOCAL",NLOCAL);
			OUT(ofs_running,"NBASE",NLOCAL);
			OUT(ofs_running,"NBANDS",NBANDS);
		}
	}

	timer::tick("Run_Frag","init",'B');
    return;
}
