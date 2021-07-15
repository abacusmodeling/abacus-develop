#include "run_md_classic.h"
#include "MD_basic.h"
#include "../src_pw/global.h"
#include "../src_io/print_info.h"
#include "../module_neighbor/sltk_atom_arrange.h"

Run_MD_CLASSIC::Run_MD_CLASSIC():grid_neigh(test_deconstructor, test_grid_driver, test_grid)
{
    pos_old1 = new double[1];
	pos_old2 = new double[1];
	pos_now = new double[1];
	pos_next = new double[1];
}

Run_MD_CLASSIC::~Run_MD_CLASSIC()
{
    delete[] pos_old1;
	delete[] pos_old2;
	delete[] pos_now;
	delete[] pos_next;
}

void Run_MD_CLASSIC::classic_md_line(void)
{
	// Setup the unitcell.
    ucell.setup_cell_classic(global_atom_card, out, ofs_running);
	DONE(ofs_running, "SETUP UNITCELL");

	Run_MD_CLASSIC run_md_classic;

	atom_arrange::search(
			SEARCH_PBC,
			ofs_running,
			run_md_classic.grid_neigh,
			ucell, 
			SEARCH_RADIUS, 
			test_atom_input,
			INPUT.test_just_neighbor);

	Print_Info PI;
    PI.setup_parameters();

	run_md_classic.md_cells_classic();

	return;
}

void Run_MD_CLASSIC::md_cells_classic(void)
{
    TITLE("Run_MD_CLASSIC", "md_cells_classic");
    timer::tick("Run_MD_CLASSIC", "md_cells_classic");

    this->md_allocate_ions();

    MD_basic mdb(INPUT.mdp, ucell);
    int mdtype = INPUT.mdp.mdtype;

    this->istep = 1;
    bool stop = false;

    while (istep <= NSTEP && !stop)
    {
        time_t estart = time(NULL);

        if (OUT_LEVEL == "ie")
        {
            cout << " -------------------------------------------" << endl;    
            cout << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
            cout << " -------------------------------------------" << endl;
            ofs_running << " -------------------------------------------" << endl;
            ofs_running << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
            ofs_running << " -------------------------------------------" << endl;
        }

		this->update_pos_classic();

		time_t eend = time(NULL);
        time_t fstart = time(NULL);

		if (mdtype == 1 || mdtype == 2)
        {
            mdb.runNVT(istep);
        }
        else if (mdtype == 0)
        {
            mdb.runNVE(istep);
        }
        else if (mdtype == -1)
        {
            stop = mdb.runFIRE(istep);
        }
        else
        {
            WARNING_QUIT("md_cells_classic", "mdtype should be -1~2!");
        }

        time_t fend = time(NULL);

		ucell.save_cartesian_position(this->pos_next);

		++istep;
    }

    timer::tick("Run_MD_CLASSIC", "md_cells_classic");
    return;
}

void Run_MD_CLASSIC::md_allocate_ions(void)
{
	pos_dim = ucell.nat * 3;

	delete[] this->pos_old1;
	delete[] this->pos_old2;
	delete[] this->pos_now;
	delete[] this->pos_next;

	this->pos_old1 = new double[pos_dim];
	this->pos_old2 = new double[pos_dim];
	this->pos_now = new double[pos_dim];
	this->pos_next = new double[pos_dim];

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
	ucell.save_cartesian_position(this->pos_now);
	return;
}