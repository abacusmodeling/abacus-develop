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

	time_t time_start = std::time(NULL);
	timer::start();

	// (1) read the input parameters.
	this->reading();

	// (2) welcome to the classic MD!
	this->classic_world();

	// (3) output information
	time_t	time_finish= std::time(NULL);
	this->print_time(time_start, time_finish);

	// (4) close all of the running logs 
	INPUT.close_log();

	return;
}

void Driver_classic::reading(void)
{
	timer::tick("Driver_classic", "reading");

	// (1) read INPUT 
	INPUT.Init( GlobalV::global_in_card );

    // (2) Print the parameters into INPUT file.
    std::stringstream ss1;
    ss1 << GlobalV::global_out_dir << GlobalV::global_in_card;
    INPUT.Print( ss1.str() );

	timer::tick("Driver_classic","reading");
	return;
}

void Driver_classic::convert(UnitCell_pseudo &ucell_c)
{
    TITLE("Driver_classic","convert");
	timer::tick("Driver_classic","convert");

    if(INPUT.atom_file!="") GlobalV::global_atom_card = INPUT.atom_file;
    GlobalV::CALCULATION = INPUT.calculation;
    GlobalV::OUT_LEVEL = INPUT.out_level;
    GlobalV::SEARCH_RADIUS = max(INPUT.search_radius,INPUT.mdp.rcut_lj+2*ANGSTROM_AU);
	GlobalV::SEARCH_PBC = INPUT.search_pbc;
    GlobalV::NSTEP = INPUT.nstep;

	GlobalV::FORCE = INPUT.force;
    GlobalV::FORCE_THR = INPUT.force_thr;

	GlobalV::STRESS = INPUT.stress;
    GlobalV::STRESS_THR = INPUT.stress_thr;
    GlobalV::PRESS1 = INPUT.press1;
    GlobalV::PRESS2 = INPUT.press2;
    GlobalV::PRESS3 = INPUT.press3;

    ucell_c.latName = INPUT.latname; 
	ucell_c.ntype = INPUT.ntype;
	ucell_c.set_vel = INPUT.set_vel;

}

void Driver_classic::classic_world(void)
{
	TITLE("Driver_classic", "classic_world");

	Run_MD_CLASSIC run_md_classic;

	this->convert(run_md_classic.ucell_c);

	if(GlobalV::CALCULATION != "md")
	{
		WARNING_QUIT("Driver_classic::classic_world","CALCULATION must be md!");
	}

    run_md_classic.classic_md_line();

	timer::finish( GlobalV::ofs_running );

	Memory::print_all( GlobalV::ofs_running ) ;

	return;
}

void Driver_classic::print_time(time_t &time_start, time_t &time_finish)
{
    // print out information before ABACUS-CMD ends
	std::cout << "\n START  Time  : " << ctime(&time_start);
	std::cout << " FINISH Time  : " << ctime(&time_finish);
	std::cout << " TOTAL  Time  : " << difftime(time_finish, time_start) << std::endl;
	std::cout << " SEE INFORMATION IN : " << GlobalV::global_out_dir << std::endl;

	GlobalV::ofs_running << "\n Start  Time  : " << ctime(&time_start);
	GlobalV::ofs_running << " Finish Time  : " << ctime(&time_finish);

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	GlobalV::ofs_running << " Total  Time  : " << hour << " h "
	    << mins << " mins "
	    << secs << " secs "<< std::endl;
}