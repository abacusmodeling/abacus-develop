#include "global.h"
#include "tools.h"
#include "cal_test0.h"
#include "../src_global/sltk_atom_arrange.h"
//xiaohui modified 2013-03-23, adding "//" before #include...
//#include "../src_develop/src_siao/trace_rho_hs.h"

int Cal_Test::Hnnz=0;

void Cal_Test::adjacent_atoms(void)
{
	TITLE("Cal_Test","adjacent_atoms");

	// (1) Find adjacent atoms for each atom.
	atom_arrange::set_sr_NL();
	atom_arrange::search( SEARCH_RADIUS );

	return;
}
// xiaohui modified 2013-03-23
/*
void Cal_Test::sparsity(void)
{
	TITLE("Cal_Test","sparsity");

	Trace_Rho_HS TRHS;
	// nonzero need to be update later.
	bool** nonzero = new bool*[NLOCAL];
	for(int col=0; col<NLOCAL; ++col)
	{
		nonzero[col] = new bool[NLOCAL-col];
		for(int row=col; row<NLOCAL; ++row)
		{
			int index = row-col;
			nonzero[col][index] = false;
			//nonzero[col][index] = true;
		}
	}
	Hnnz=TRHS.cal_Hnnz(nonzero);
	OUT(ofs_running,"Hnnz",Hnnz);
	cout << " Hnnz = " << Hnnz << endl;

	for(int col=0; col<NLOCAL; ++col)
	{
		delete[] nonzero[col];
	}
	delete[] nonzero;
	
	return;
}
*/
