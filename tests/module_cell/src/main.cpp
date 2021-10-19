//#include "timer.h"
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "cell_unittest.h"

using namespace std;

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Hello, this is the cell module of ABACUS." << std::endl;

	calculate();

    return 0;
}

void calculate()
{
#ifdef __LCAO
	test_cell_orb test;
#else
	test_cell test;
#endif

	test.count_ntype();
	test.set_parameters();
	test.setup_cell();

#ifdef __LCAO
	test.setup_kpt();
	test.set_ekcut();
	test.set_orbs(test.ucell.lat0);
	test.prep_neighbour();
	test.alloc_ST();
	test.build_ST_new('S');
	test.build_ST_new('T');
	test.print_ST();
#endif
    return;
}
