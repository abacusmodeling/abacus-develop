#include "driver_classic.h"
#include "run_md_classic.h"
#include "../input.h"
#include "../module_base/tool_title.h"
#include "../module_base/timer.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"

Driver_classic::Driver_classic(){}
Driver_classic::~Driver_classic(){}

void Driver_classic::init()
{
	TITLE("Driver_classic", "init");

	// (1) read the input parameters.
	this->reading();

	// (2) welcome to the classic MD!
	this->classic_world();

	// (3) close all of the running logs 
	INPUT.close_log();

	return;
}

void Driver_classic::reading(void)
{
	timer::tick("Driver_classic", "reading");

	// (1) read INPUT 
	INPUT.Init( global_in_card );

	// (2) copy the variables from INPUT to each class
	this->convert();

    // (3) Print the parameters into INPUT file.
    stringstream ss1;
    ss1 << global_out_dir << global_in_card;
    INPUT.Print( ss1.str() );

	timer::tick("Driver_classic","reading");
	return;
}

void Driver_classic::convert(void)
{
    TITLE("Driver_classic","convert");
	timer::tick("Driver_classic","convert");

    if(INPUT.atom_file!="") global_atom_card = INPUT.atom_file;
    CALCULATION = INPUT.calculation;
    OUT_LEVEL = INPUT.out_level;
    SEARCH_RADIUS = INPUT.search_radius;
	SEARCH_PBC = INPUT.search_pbc;
    NSTEP = INPUT.nstep;

    ucell.latName = INPUT.latname; 
	ucell.ntype = INPUT.ntype;

}

void Driver_classic::classic_world(void)
{
	TITLE("Driver_classic", "classic_world");

	Run_MD_CLASSIC run_md_classic;
    run_md_classic.classic_md_line();

	timer::finish( ofs_running );

	Memory::print_all( ofs_running ) ;

	return;
}

