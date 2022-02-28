#include "run_lcao.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "module_neighbor/sltk_atom_arrange.h"
#include "src_lcao/LOOP_cell.h"
#include "src_io/print_info.h"
#include "src_lcao/run_md_lcao.h"

Run_lcao::Run_lcao(){}
Run_lcao::~Run_lcao(){}


void Run_lcao::lcao_line(ModuleEnSover::En_Solver *p_ensolver)
{
    ModuleBase::TITLE("Run_lcao","lcao_line");
    ModuleBase::timer::tick("Run_lcao", "lcao_line");
    //------------------------------------init Cell----------------------------
    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
#ifdef __LCAO
    GlobalC::ucell.setup_cell( GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::global_atom_card, GlobalV::ofs_running);
#else
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalV::global_atom_card, GlobalV::ofs_running);
#endif
	if(INPUT.test_just_neighbor)
	{
		//test_search_neighbor();
		GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
			GlobalV::ofs_running,
			GlobalV::OUT_LEVEL,
			GlobalC::ORB.get_rcutmax_Phi(),
			GlobalC::ucell.infoNL.get_rcutmax_Beta(),
			GlobalV::GAMMA_ONLY_LOCAL);

		atom_arrange::search(
			GlobalV::SEARCH_PBC,
			GlobalV::ofs_running,
			GlobalC::GridD,
			GlobalC::ucell,
			GlobalV::SEARCH_RADIUS,
			GlobalV::test_atom_input,
			INPUT.test_just_neighbor);
    }
    // the symmetry of a variety of systems.
    if (GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        ModuleBase::QUIT();
    }
    //------------------------------------\init Cell----------------------------


    //------------------------------------------------------------
    //---------------------Init En_Solver-------------------------
    p_ensolver->Init(INPUT, GlobalC::ucell);
    //------------------------------------------------------------

    //------------------------------------init Basis_lcao----------------------------
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    ORB_control orb_con(
        GlobalV::GAMMA_ONLY_LOCAL,
        GlobalV::NLOCAL, GlobalV::NBANDS,
        GlobalV::NSPIN, GlobalV::DSIZE,
        GlobalV::NB2D, GlobalV::DCOLOR,
        GlobalV::DRANK, GlobalV::MY_RANK,
        GlobalV::CALCULATION, GlobalV::KS_SOLVER);

    orb_con.read_orb_first(
		GlobalV::ofs_running,
		GlobalC::ORB,
		GlobalC::ucell.ntype,
		GlobalC::ucell.lmax,
		INPUT.lcao_ecut,
		INPUT.lcao_dk,
		INPUT.lcao_dr,
		INPUT.lcao_rmax,
		GlobalV::out_descriptor,
		INPUT.out_r_matrix,
		GlobalV::FORCE,
		GlobalV::MY_RANK);
		
	GlobalC::ucell.infoNL.setupNonlocal(
		GlobalC::ucell.ntype,
		GlobalC::ucell.atoms,
		GlobalV::ofs_running,
		GlobalC::ORB
	);

#ifdef __MPI   
	orb_con.set_orb_tables(
		GlobalV::ofs_running,
		GlobalC::UOT,
		GlobalC::ORB,
		GlobalC::ucell.lat0,
		GlobalV::out_descriptor,
		Exx_Abfs::Lmax,
		GlobalC::ucell.infoNL.nprojmax,
		GlobalC::ucell.infoNL.nproj,
        GlobalC::ucell.infoNL.Beta);
#endif

    orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);
    //------------------------------------\init Basis_lcao----------------------------

    if (GlobalV::CALCULATION == "md")
	{
		Run_MD_LCAO run_md_lcao(orb_con.ParaV);
		run_md_lcao.opt_cell(orb_con, p_ensolver);
	}
	else // cell relaxations
	{
        LOOP_cell lc(orb_con.ParaV);
        //keep wfc_gamma or wfc_k remaining
        lc.opt_cell(orb_con, p_ensolver);
	}

	ModuleBase::timer::tick("Run_lcao","lcao_line");
    return;
}
