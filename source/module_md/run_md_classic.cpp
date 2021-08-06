#include "run_md_classic.h"
#include "MD_basic.h"
#include "LJ_potential.h"
#include "../input.h"
#include "../module_base/global_variable.h"
#include "../src_io/print_info.h"
#include "../module_neighbor/sltk_atom_arrange.h"

Run_MD_CLASSIC::Run_MD_CLASSIC():grid_neigh(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid)
{
    pos_old1 = new double[1];
	pos_old2 = new double[1];
	pos_now = new double[1];
	pos_next = new double[1];
	force=new Vector3<double>[1];
	stress.create(3,3);
}

Run_MD_CLASSIC::~Run_MD_CLASSIC()
{
    delete[] pos_old1;
	delete[] pos_old2;
	delete[] pos_now;
	delete[] pos_next;
	delete[] force;
}

void Run_MD_CLASSIC::classic_md_line(void)
{
	// Setup the unitcell.
    ucell_c.setup_cell_classic(GlobalV::global_atom_card, GlobalV::ofs_running, GlobalV::ofs_warning);
	DONE(GlobalV::ofs_running, "SETUP UNITCELL");

	//Print_Info PI;
    //PI.setup_parameters();

	this->md_cells_classic();

	return;
}

void Run_MD_CLASSIC::md_cells_classic(void)
{
    TITLE("Run_MD_CLASSIC", "md_cells_classic");
    timer::tick("Run_MD_CLASSIC", "md_cells_classic");

    this->md_allocate_ions();

    MD_basic mdb(INPUT.mdp, ucell_c);
    int mdtype = INPUT.mdp.mdtype;

    this->istep = 1;
    bool stop = false;

    while (istep <= GlobalV::NSTEP && !stop)
    {
		time_t fstart = time(NULL);

        if (GlobalV::OUT_LEVEL == "ie")
        {
            cout << " -------------------------------------------" << endl;    
            cout << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
            cout << " -------------------------------------------" << endl;
            GlobalV::ofs_running << " -------------------------------------------" << endl;
            GlobalV::ofs_running << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
            GlobalV::ofs_running << " -------------------------------------------" << endl;
        }

		double potential;

		if(this->ucell_c.judge_big_cell())
		{
			
		
		}
		else
		{
			atom_arrange::search(
				GlobalV::SEARCH_PBC,
				GlobalV::ofs_running,
				this->grid_neigh,
				this->ucell_c, 
				GlobalV::SEARCH_RADIUS, 
				GlobalV::test_atom_input,
				INPUT.test_just_neighbor);

			LJ_potential LJ_CMD;

			potential = LJ_CMD.Lennard_Jones(
							this->ucell_c,
							this->grid_neigh,
							this->force,
							this->stress);
		}

		this->update_pos_classic();

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
            WARNING_QUIT("md_cells_classic", "mdtype should be -1~2!");
        }

        time_t fend = time(NULL);

		ucell_c.save_cartesian_position(this->pos_next);

		++istep;
    }

    timer::tick("Run_MD_CLASSIC", "md_cells_classic");
    return;
}

void Run_MD_CLASSIC::md_allocate_ions(void)
{
	pos_dim = ucell_c.nat * 3;

	delete[] this->pos_old1;
	delete[] this->pos_old2;
	delete[] this->pos_now;
	delete[] this->pos_next;
	delete[] this->force;

	this->pos_old1 = new double[pos_dim];
	this->pos_old2 = new double[pos_dim];
	this->pos_now = new double[pos_dim];
	this->pos_next = new double[pos_dim];
	this->force = new Vector3<double>[ucell_c.nat];

	ZEROS(pos_old1, pos_dim);
	ZEROS(pos_old2, pos_dim);
	ZEROS(pos_now, pos_dim);
	ZEROS(pos_next, pos_dim);
}

void Run_MD_CLASSIC::update_pos_classic(void)
{
	for(int i=0;i<pos_dim;i++)
	{
		this->pos_old2[i] = this->pos_old1[i];
		this->pos_old1[i] = this->pos_now[i];
	}
	this->ucell_c.save_cartesian_position(this->pos_now);
	return;
}