#include "unitcell_pseudo.h"

#ifdef __LCAO
#include "../module_orbital/ORB_unittest.h"
#include "../module_orbital/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30

class test_cell_orb : public test_orb
#else
class test_cell
#endif
{

public:
#ifdef __LCAO
	test_cell_orb();
	~test_cell_orb();
#else
	test_cell();
	~test_cell();
#endif

	std::ofstream ofs_warning;
	UnitCell_pseudo ucell;
	output out;

#ifndef __LCAO
	std::ofstream ofs_running;
	int ntype;
#endif

	void set_parameters();//set some global variables
	void setup_cell();

#ifndef __LCAO
	void count_ntype(); //from STRU, count types of elements
#endif
};

