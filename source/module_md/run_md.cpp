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
#include "../src_pw/symmetry.h"
#include "../module_neighbor/sltk_atom_arrange.h"

Run_md::Run_md(){}
Run_md::~Run_md(){}


void Run_md::md_line(void)
{
	TITLE("Run_md","md_line");
	timer::tick("Run_md","md_line");

		Run_md::classic_md_line();

	timer::tick("Run_md","md_line");
    return;
}

void Run_md::classic_md_line(void)
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