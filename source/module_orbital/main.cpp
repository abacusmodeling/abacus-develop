//#include "timer.h"
#include <ctime>
#include "ORB_control.h"
#include "ORB_unittest.h"
#include "../module_base/global_variable.h"

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Hello, this is the ORB module of ABACUS." << std::endl;

    calculate();

    return 0;
}


void calculate()
{
	GlobalV::ofs_running.open("log.txt");
	GlobalV::ofs_warning.open("warning.txt");
	GlobalV::BASIS_TYPE = "lcao";

	test_orb test;

	test.count_ntype();
	test.set_files();
/*
	ORB_control ooo;
	ooo.Read_Orbitals(
	GlobalV::ofs_running,
	test.ntype,
	test.lmax,
	test.out_descriptor,
	test.out_r_matrix,
	test.force_flag,
	test.my_rank);

*/
//	ooo.set_orb_tables();


	std::cout << "--------------------" << std::endl;
	std::cout << " Have a great day! " << std::endl;
	std::cout << "--------------------" << std::endl;

	GlobalV::ofs_running.close();
	GlobalV::ofs_warning.close();

    return;
}
