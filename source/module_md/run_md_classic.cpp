#include "run_md_classic.h"
#include "MD_basic.h"
#include "LJ_potential.h"
#include "DP_potential.h"
#include "cmd_neighbor.h"
#include "../input.h"
#include "../src_io/print_info.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_neighbor/sltk_grid_driver.h"

Run_MD_CLASSIC::Run_MD_CLASSIC()
{
	force = new ModuleBase::Vector3<double>[1];
	stress.create(3,3);
}

Run_MD_CLASSIC::~Run_MD_CLASSIC()
{
	delete[] force;
}

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

    this->md_allocate_ions();

    MD_basic mdb(INPUT.mdp, ucell_c);
    int mdtype = INPUT.mdp.mdtype;

    int istep = 1;
    bool stop = false;
	double potential;

    while (istep <= GlobalV::NSTEP && !stop)
    {
		time_t fstart = time(NULL);

        Print_Info::print_screen(0, 0, istep);

		this->md_force_stress(potential);

		if (mdtype == 1 || mdtype == 2)
        {
            mdb.runNVT(istep, potential, force, stress);
        }
        else if (mdtype == 0)
        {
            mdb.runNVE(istep, potential, force, stress);
        }
        else if (mdtype == -1)
        {
            stop = mdb.runFIRE(istep, potential, force, stress);
        }
        else
        {
            ModuleBase::WARNING_QUIT("md_cells_classic", "mdtype should be -1~2!");
        }

        time_t fend = time(NULL);

		++istep;
    }

	GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << potential*ModuleBase::Hartree_to_eV << " eV" << std::endl; 
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    ModuleBase::timer::tick("Run_MD_CLASSIC", "md_cells_classic");
    return;
}

void Run_MD_CLASSIC::md_force_stress(double &potential)
{
	ModuleBase::TITLE("Run_MD_CLASSIC", "md_force_stress");
    ModuleBase::timer::tick("Run_MD_CLASSIC", "md_force_stress");

	if(INPUT.mdp.md_potential == "LJ")
	{
		bool which_method = this->ucell_c.judge_big_cell();
		if(which_method)
		{
			CMD_neighbor cmd_neigh;
			cmd_neigh.neighbor(this->ucell_c);

			potential = LJ_potential::Lennard_Jones(
								this->ucell_c,
								cmd_neigh,
								this->force,
								this->stress);
		}
		else
		{
			Grid_Driver grid_neigh(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid);
			atom_arrange::search(
					GlobalV::SEARCH_PBC,
					GlobalV::ofs_running,
					grid_neigh,
					this->ucell_c, 
					GlobalV::SEARCH_RADIUS, 
					GlobalV::test_atom_input,
					INPUT.test_just_neighbor);

			potential = LJ_potential::Lennard_Jones(
								this->ucell_c,
								grid_neigh,
								this->force,
								this->stress);
		}
	}
	else if(INPUT.mdp.md_potential == "DP")
	{
		DP_potential::DP_pot(ucell_c, potential, force, stress);
	}
	else if(INPUT.mdp.md_potential == "FP")
	{
		ModuleBase::WARNING_QUIT("md_force_stress", "FPMD is only available in integrated program or PW module ！");
	}
	else
	{
		ModuleBase::WARNING_QUIT("md_force_stress", "Unsupported MD potential ！");
	}

	ModuleBase::GlobalFunc::NEW_PART("   TOTAL-FORCE (eV/Angstrom)");
	GlobalV::ofs_running << std::endl;
	GlobalV::ofs_running << " atom    x              y              z" << std::endl;
	const double fac = ModuleBase::Hartree_to_eV*ModuleBase::ANGSTROM_AU;
	int iat = 0;
	for(int it=0; it<ucell_c.ntype; ++it)
	{
		for(int ia=0; ia<ucell_c.atoms[it].na; ++ia)
		{
			std::stringstream ss;
			ss << ucell_c.atoms[it].label << ia+1;
			GlobalV::ofs_running << " " << std::left << std::setw(8) << ss.str()
						  		<< std::setw(15) << std::setiosflags(ios::fixed) << std::setprecision(6) << force[iat].x*fac
						  		<< std::setw(15) << std::setiosflags(ios::fixed) << std::setprecision(6) << force[iat].y*fac
						  		<< std::setw(15) << std::setiosflags(ios::fixed) << std::setprecision(6) << force[iat].z*fac
						  		<< std::endl;
			iat++;
		}
	}

	ModuleBase::timer::tick("Run_MD_CLASSIC", "md_force_stress");
}


void Run_MD_CLASSIC::md_allocate_ions(void)
{
	int pos_dim = ucell_c.nat * 3;

	delete[] this->force;
	this->force = new ModuleBase::Vector3<double>[ucell_c.nat];
}
