#include "../../../source/module_cell/unitcell_pseudo.h"
//#include "../../../source/src_parallel/parallel_global.h"
#include "../../../source/module_base/global_variable.h"

#ifdef __LCAO
#include "../../module_orb/src/ORB_unittest.h"
#include "../../../source/module_orbital/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#include "../../../source/module_neighbor/sltk_grid_driver.h"
#include "../../../source/module_neighbor/sltk_atom_arrange.h"
#include "klist.h"
#include <sstream>

namespace Test_Cell
{
	extern Grid_Driver GridD;
}

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

#ifdef __LCAO
	int GAMMA_ONLY_LOCAL = 1;
	Test_Cell::K_Vectors kv;
	double **Sgamma;
	double **Tgamma;
	std::complex<double> ***Sk;
	std::complex<double> ***Tk;
#else
	std::ofstream ofs_running;
	int ntype;
#endif

	void set_parameters();//set some global variables
	void setup_cell();

#ifdef __LCAO
	void setup_kpt();
	void prep_neighbour();
	void alloc_ST();
	void build_ST_new(const char& dtype);
	void print_ST();
#else
	void count_ntype(); //from STRU, count types of elements
#endif
};

