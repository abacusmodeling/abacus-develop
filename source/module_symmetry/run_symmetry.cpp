#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "symmetry.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_global.h"
using namespace std;

void calculate();

int main(int argc, char **argv)
{

	cout << "Hello, this is the 'symmetry' module of ABACUS." << endl;

	cout << "The module does symmetry analysis for an input geometry." << endl;
	Parallel_Global::read_mpi_parameters(argc,argv);
	//cout << "Right now, the module is still empty, soon we will have more tests." << endl;

    calculate();

    return 0;
}


void calculate()
{
	//ofstream ofs("log.txt");
	ofstream ofs_running("log.txt");
	ofstream ofs("useless.txt");
	ofstream ofs_warning("warning.txt");
	ifstream ifs("INPUT");
	UnitCell_pseudo ucell;
	Symmetry symm;
	ifs >> ucell.ntype;
	ifs.close();
	output out;
	ucell.setup_cell_classic(
	"STRU", 
	out, 
	ofs,
	ofs_running,
	ofs_warning)
	symm.analy_sys(ucell, out, ofs_running);
	ofs_running.close();
//	ooo.set_orb_tables();

	//ofs.close();

	cout << "--------------------" << endl;
	cout << " Have a great day! " << endl;
	cout << "--------------------" << endl;


    return;
}
