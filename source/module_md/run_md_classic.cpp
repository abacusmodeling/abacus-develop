#include "run_md_classic.h"
#include "MD_basic.h"
#include "LJ_potential.h"
#include "DP_potential.h"
#include "cmd_neighbor.h"
#include "../src_io/print_info.h"
#include "../input.h"
#include "NVE.h"

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

	// if(INPUT.mdp.mdtype==0)
	// {
	// 	NVE verlet;
	// }

	bool stop = false;

	NVE verlet(INPUT.mdp, ucell_c);
	verlet.setup();

	return;

    // while (istep <= GlobalV::NSTEP && !stop)
    // {
	// 	time_t fstart = time(NULL);

    //     Print_Info::print_screen(0, 0, istep);

	// 	this->md_force_stress(potential);

	// 	if (mdtype == 1 || mdtype == 2)
    //     {
    //         mdb.runNVT(istep, potential, force, stress);
    //     }
    //     else if (mdtype == 0)
    //     {
    //         mdb.runNVE(istep, potential, force, stress);
    //     }
    //     else if (mdtype == -1)
    //     {
    //         stop = mdb.runFIRE(istep, potential, force, stress);
    //     }
    //     else
    //     {
    //         ModuleBase::WARNING_QUIT("md_cells_classic", "mdtype should be -1~2!");
    //     }

    //     time_t fend = time(NULL);

	// 	++istep;
    // }

	// GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    // GlobalV::ofs_running << std::setprecision(16);
    // GlobalV::ofs_running << " !FINAL_ETOT_IS " << potential*ModuleBase::Hartree_to_eV << " eV" << std::endl; 
    // GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    ModuleBase::timer::tick("Run_MD_CLASSIC", "md_cells_classic");
    return;
}