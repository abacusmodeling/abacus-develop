#include "run_pw.h"
#include "src_pw/global.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "src_io/optical.h"
#include "src_io/numerical_basis.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/print_info.h"
#include "src_ions/ions.h"
#include "src_pw/run_md_pw.h"

Run_pw::Run_pw(){}
Run_pw::~Run_pw(){}

void Run_pw::plane_wave_line(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Run_pw","plane_wave_line");
	ModuleBase::timer::tick("Run_pw","plane_wave_line");

    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
#ifdef __LCAO
    GlobalC::ucell.setup_cell( GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#else
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#endif

    // mohan add 2010-10-10, just to test the symmetry of a variety
    // of systems.
    if(GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        ModuleBase::QUIT();
    }

    //------------------------------------------------------------
    //---------------------Init ESolver-------------------------
    //------------------------------------------------------------
    p_esolver->Init(INPUT, GlobalC::ucell);

    

    if(GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "sto-md")
    {
        Run_MD_PW run_md_pw;
        run_md_pw.md_ions_pw(p_esolver);
    }
    else
    {
        Ions ions;
        ions.opt_ions_pw(p_esolver);
    }

    // cout<<"cpws SUCCESS"<<endl;


	if(Optical::opt_epsilon2 && (GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw"))
	{
		Optical opt;
		opt.cal_epsilon2(GlobalV::NBANDS);
	}

    p_esolver->postprocess();

	ModuleBase::timer::tick("Run_pw","plane_wave_line");
    return;
}
