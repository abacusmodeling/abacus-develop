#include "run_md_classic.h"
#include "MD_basic.h"
#include "NVE.h"
#include "MSST.h"
#include "NVT_ADS.h"
#include "../input.h"
#include "../src_io/print_info.h"

Run_MD_CLASSIC::Run_MD_CLASSIC(){}

Run_MD_CLASSIC::~Run_MD_CLASSIC(){}

void Run_MD_CLASSIC::classic_md_line(void)
{
	ModuleBase::TITLE("Run_MD_CLASSIC", "classic_md_line");
    ModuleBase::timer::tick("Run_MD_CLASSIC", "classic_md_line");

	// Setup the unitcell.
#ifdef __LCAO
	ucell_c.setup_cell_classic(GlobalC::ORB, GlobalV::global_atom_card, GlobalV::ofs_running, GlobalV::ofs_warning);
#else
    ucell_c.setup_cell_classic(GlobalV::global_atom_card, GlobalV::ofs_running, GlobalV::ofs_warning);
#endif
	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // determine the mdtype
    Verlet *verlet;
    if(INPUT.mdp.mdtype==0)
    {
        verlet = new NVE(INPUT.mdp, ucell_c); 
    }
    else if(INPUT.mdp.mdtype==1)
    {
        verlet = new NVT_ADS(INPUT.mdp, ucell_c);
    }
    else if(INPUT.mdp.mdtype==4)
    {
        verlet = new MSST(INPUT.mdp, ucell_c); 
    }

    bool stop = false;

    // md cycle
    while (verlet->step_ <= (GlobalV::NSTEP + verlet->step_rst_) && !stop)
    {
        if(verlet->step_ == verlet->step_rst_)
        {
            verlet->setup();
        }
        else
        {
            verlet->first_half();

            // update force and virial due to the update of atom positions
            MD_func::force_virial(verlet->mdp, verlet->ucell, verlet->potential, verlet->force, verlet->virial);

            verlet->second_half();

            MD_func::kinetic_stress(verlet->ucell, verlet->vel, verlet->allmass, verlet->kinetic, verlet->stress);

            verlet->stress +=  verlet->virial;
        }

        Print_Info::print_screen(0, 0, verlet->step_);
        verlet->outputMD();

        if((verlet->step_ - verlet->step_rst_) % verlet->mdp.recordFreq == 0)
        {
            verlet->ucell.update_vel(verlet->vel);
            std::stringstream file;
            file << GlobalV::global_out_dir << "STRU_MD_" << verlet->step_;
#ifdef __LCAO
            verlet->ucell.print_stru_file(GlobalC::ORB, file.str(), 1, 1);
#else
            verlet->ucell.print_stru_file(file.str(), 1, 1);
#endif
            MD_func::MDdump(verlet->step_, verlet->ucell.nat, verlet->virial, verlet->force);
            verlet->write_restart();
        }

        verlet->step_++;
    }

	GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << verlet->potential*ModuleBase::Hartree_to_eV << " eV" << std::endl; 
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    delete verlet;
    ModuleBase::timer::tick("Run_MD_CLASSIC", "md_cells_classic");
    return;
}