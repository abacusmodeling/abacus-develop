#include "run_md.h"
#include "run_md_pw.h"
#ifdef __LCAO
#include "run_md_lcao.h"
#endif
#include "run_md_classic.h"
#include "../src_pw/global.h"
#include "../input.h"
#include "../src_io/cal_test.h"
#include "../src_io/print_info.h"
#include "../module_symmetry/symmetry.h"

Run_md::Run_md(){}
Run_md::~Run_md(){}


void Run_md::md_line(void)
{
	TITLE("Run_md","md_line");
	timer::tick("Run_md","md_line");

	if(INPUT.mdp.md_potential==0)
	{
		Run_md::ai_md_line();
	}
	else
	{
		Run_md::classic_md_line();
	}

	timer::tick("Run_md","md_line");
    return;
}


void Run_md::ai_md_line(void)
{
    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalC::out, GlobalV::global_atom_card, GlobalV::ofs_running);

	// setup GlobalV::NBANDS 
	// Yu Liu add 2021-07-03
	GlobalC::CHR.cal_nelec();

	// mohan add 2010-09-06
	// Yu Liu move here 2021-06-27
	// because the number of element type
	// will easily be ignored, so here
	// I warn the user again for each type.
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		GlobalC::xcf.which_dft(GlobalC::ucell.atoms[it].dft);
	}

    //GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir , GlobalV::global_atom_card , GlobalV::ofs_running, GlobalV::NLOCAL, GlobalV::NBANDS);
    DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (Symmetry::symm_flag)
    {
        symm.analy_sys(GlobalC::ucell, GlobalC::out, GlobalV::ofs_running);
        DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set(symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec );
    DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info PI;
    PI.setup_parameters();

#ifdef __LCAO
    // * reading the localized orbitals/projectors 
	// * construct the interpolation tables.
	if(GlobalV::BASIS_TYPE=="lcao")
	{
		LOWF.orb_con.set_orb_tables(
		GlobalV::ofs_running,
		UOT, 
		ORB,
		GlobalC::ucell.ntype,
		GlobalC::ucell.lmax,
		INPUT.lcao_ecut,
		INPUT.lcao_dk,
		INPUT.lcao_dr,
		INPUT.lcao_rmax, 
		GlobalC::ucell.lat0, 
		INPUT.out_descriptor,
		INPUT.out_r_matrix,
		Exx_Abfs::Lmax,
		GlobalV::FORCE,
		GlobalV::MY_RANK);

		// * allocate H and S matrices according to computational resources
		// * set the 'trace' between local H/S and global H/S
		LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, ParaO);
	}
#endif

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, GlobalC::ucell, GlobalC::kv);
    DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << endl;
    cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << endl;

    // the symmetry of a variety of systems.
    if(GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        QUIT();
    }

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
        GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04

	// Inititlize the charge density.
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    DONE(GlobalV::ofs_running,"INIT CHARGE");

	// Initializee the potential.
    GlobalC::pot.allocate(GlobalC::pw.nrxx);
    DONE(GlobalV::ofs_running,"INIT POTENTIAL");

	if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		Run_MD_PW run_md_pw;
        run_md_pw.md_cells_pw();
	}
#ifdef __LCAO
	else if(GlobalV::BASIS_TYPE=="lcao")
	{
		Run_MD_LCAO run_md_lcao;
		run_md_lcao.opt_cell();
	}
#endif

    return;
}

void Run_md::classic_md_line(void)
{
	// Setup the unitcell.
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalC::out, GlobalV::global_atom_card, GlobalV::ofs_running);
	DONE(GlobalV::ofs_running, "SETUP UNITCELL");

	Print_Info PI;
    PI.setup_parameters();

	Run_MD_CLASSIC run_md_classic;
	run_md_classic.md_cells_classic();

	return;
}