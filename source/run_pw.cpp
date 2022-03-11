#include "run_pw.h"
#include "src_pw/global.h"
#include "src_pw/energy.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "src_io/numerical_basis.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/print_info.h"
#include "module_symmetry/symmetry.h"
#include "src_ions/Cell_PW.h"
#include "src_pw/run_md_pw.h"

Run_pw::Run_pw(){}
Run_pw::~Run_pw(){}

void Run_pw::plane_wave_line(void)
{
    ModuleBase::TITLE("Run_pw","plane_wave_line");
	ModuleBase::timer::tick("Run_pw","plane_wave_line");

    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
#ifdef __LCAO
    GlobalC::ucell.setup_cell( GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::global_atom_card, GlobalV::ofs_running);
#else
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalV::global_atom_card, GlobalV::ofs_running);
#endif
    //GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir , GlobalV::global_atom_card , GlobalV::ofs_running, GlobalV::NLOCAL, GlobalV::NBANDS);

    // setup GlobalV::NBANDS 
	// Yu Liu add 2021-07-03
	GlobalC::CHR.cal_nelec();

	if(GlobalC::ucell.atoms[0].xc_func=="HSE"||GlobalC::ucell.atoms[0].xc_func=="PBE0")
	{
		XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
	}

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set( GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(GlobalC::ucell, GlobalC::kv);

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, GlobalC::ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << std::endl;

    // mohan add 2010-10-10, just to test the symmetry of a variety
    // of systems.
    if(GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        ModuleBase::QUIT();
    }

    // mohan add 2010-09-13
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
        GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04
        

    // Calculate Structure factor
    GlobalC::pw.setup_structure_factor();
    // cout<<"after pgrid init nrxx = "<<GlobalC::pw.nrxx<<endl;
    
//----------------------------------------------------------
// 1 read in initial data:
//   a lattice structure:atom_species,atom_positions,lattice vector
//   b k_points
//   c pseudopotential
// 2 setup planeware basis, FFT,structure factor, ...
// 3 initialize local and nonlocal pseudopotential in G_space
// 4 initialize charge desity and warefunctios in G_space
//----------------------------------------------------------

    //=====================================
    // init charge/potential/wave functions
    //=====================================
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    GlobalC::pot.allocate(GlobalC::pw.nrxx);

    if(GlobalV::CALCULATION == "md")
    {
        Run_MD_PW run_md_pw;
        run_md_pw.md_cells_pw();
    }
    else
    {
        Cell_PW cpws;
        cpws.opt_cells_pw();
    }

    // cout<<"cpws SUCCESS"<<endl;


    // caoyu add 2020-11-24, mohan updat 2021-01-03
    if(GlobalV::BASIS_TYPE=="pw" && GlobalV::out_descriptor==1)
    {
        Numerical_Descriptor nc;
        nc.output_descriptor(GlobalC::wf.evc, INPUT.lmax_descriptor);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"GENERATE DESCRIPTOR FOR DEEPKS");
    }


    if(GlobalV::BASIS_TYPE=="pw" && winput::out_spillage) //xiaohui add 2013-09-01
    {
        //std::cout << "\n Output Spillage Information : " << std::endl;
        // calculate spillage value.
#ifdef __LCAO
        if ( winput::out_spillage == 3)
        {
            GlobalV::BASIS_TYPE="pw"; 
            std::cout << " NLOCAL = " << GlobalV::NLOCAL << std::endl;

            for (int ik=0; ik<GlobalC::kv.nks; ik++)
            {
                GlobalC::wf.wanf2[ik].create(GlobalV::NLOCAL, GlobalC::wf.npwx);
				if(GlobalV::BASIS_TYPE=="pw")
                {
					std::cout << " ik=" << ik + 1 << std::endl;

                    GlobalV::BASIS_TYPE="lcao_in_pw";
					GlobalC::wf.LCAO_in_pw_k(ik, GlobalC::wf.wanf2[ik]);
                    GlobalV::BASIS_TYPE="pw";
                }
            }

            //Spillage sp;
            //sp.get_both(GlobalV::NBANDS, GlobalV::NLOCAL, GlobalC::wf.wanf2, GlobalC::wf.evc);
        }
#endif

        // output overlap
        if ( winput::out_spillage <= 2 )
        {
            Numerical_Basis numerical_basis;
            numerical_basis.output_overlap(GlobalC::wf.evc);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"BASIS OVERLAP (Q and S) GENERATION.");
        }
    }


	if(Optical::opt_epsilon2 && (GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw"))
	{
		Optical opt;
		opt.cal_epsilon2(GlobalV::NBANDS);
	}

	// compute density of states
	GlobalC::en.perform_dos_pw();

	ModuleBase::timer::tick("Run_pw","plane_wave_line");
    return;
}
