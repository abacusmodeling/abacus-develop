#ifndef PARALLEL_ATOMS_H
#define PARALLEL_ATOMS_H

#include "../src_pw/tools.h"

class Parallel_Atoms
{
	public:

	Parallel_Atoms();
	~Parallel_Atoms();

	bool* keep_this_atom;
	int nat; // how many atoms on this processor.
	int nlocal; // how many orbitals on this processor.
	void cut_atoms(void);	
	void set_trace(int *trace_loc_row, int *trace_loc_col, int &nrow, int &ncol);

	private:


};

#endif
