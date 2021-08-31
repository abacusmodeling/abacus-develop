//#include "timer.h"
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../src_parallel/parallel_global.h"
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
	test.set_ekcut();
	test.read_files();
#endif
    return;
}
