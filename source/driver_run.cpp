#include "driver.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/winput.h"
#include "module_neighbor/sltk_atom_arrange.h"
#include "src_lcao/LOOP_ions.h"
#include "src_io/print_info.h"
#include "src_lcao/run_md_lcao.h"
#include "src_pw/run_md_pw.h"

void Driver::driver_run()
{
    ModuleBase::TITLE("Driver", "driver_line");
    ModuleBase::timer::tick("Driver", "driver_line");

    //-----------------------init Cell--------------------------
    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters

    // 1. Initialzie type of Esolver
    ModuleESolver::ESolver *p_esolver = nullptr;
    ModuleESolver::init_esolver(p_esolver);

    // 2. Setup cell and atom information
#ifdef __LCAO
    GlobalC::ucell.setup_cell(GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#else
    if(GlobalV::BASIS_TYPE == "lcao_in_pw" || GlobalV::BASIS_TYPE == "lcao")
    {
        ModuleBase::WARNING_QUIT("driver","to use LCAO basis, compile with __LCAO");
    }
    GlobalC::ucell.setup_cell(GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#endif

    // 3. For these two types of calculations
    // nothing else need to be initialized
    if(GlobalV::CALCULATION == "test_neighbour" || GlobalV::CALCULATION == "test_memory")
    {
        p_esolver->Run(0, GlobalC::ucell);
        ModuleBase::QUIT();
    }

    // 4. Initialize Esolver
    p_esolver->Init(INPUT, GlobalC::ucell);

    if(GlobalV::BASIS_TYPE=="lcao" && GlobalV::CALCULATION=="get_S")
    {
        p_esolver->Run(0, GlobalC::ucell);
        ModuleBase::timer::tick("Driver_run", "driver_line");
        return;
    }
    //------------------------------------------------------------

    // This part onward needs to be refactored.
    //---------------------------MD/Relax------------------
    if(GlobalV::BASIS_TYPE=="lcao")
    {
        if (GlobalV::CALCULATION == "md")
        {
            Run_MD_LCAO run_md_lcao;
            run_md_lcao.opt_ions(p_esolver);
        }
        else // cell relaxations
        {
            LOOP_ions ions; 
            ions.opt_ions(p_esolver);
        }
    }
    else
    {
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

        if(Optical::opt_epsilon2)
        {
            Optical opt;
            opt.cal_epsilon2(GlobalV::NBANDS);            
        }

        p_esolver->postprocess();
    }
    //---------------------------MD/Relax------------------

    ModuleESolver::clean_esolver(p_esolver);

    ModuleBase::timer::tick("Driver", "driver_line");
    return;
}
