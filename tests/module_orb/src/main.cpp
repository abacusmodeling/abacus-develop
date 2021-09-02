//#include "timer.h"
#include <fstream>
#include <ctime>
#include "ORB_unittest.h"
#include "../../../source/module_base/global_variable.h"
#include "../../../source/module_base/global_file.h"

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Hello, this is the ORB module of ABACUS." << std::endl;

    calculate();

    return 0;
}

void calculate()
{
	GlobalV::BASIS_TYPE = "lcao";

	test_orb test;

	test.ofs_running.open("log.txt");
	test.count_ntype();
	test.set_files();
	test.set_ekcut();
	test.set_orbs(test.lat0);

	std::cout << "--------------------" << std::endl;
	std::cout << " Have a great day! " << std::endl;
	std::cout << "--------------------" << std::endl;

    return;
}
