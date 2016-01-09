#ifndef PARALLEL_ORBITALS_H
#define PARALLEL_ORBITALS_H

#include "../src_pw/tools.h"
#include "../src_external/src_pdiag/pdiag_double.h"

class Parallel_Orbitals : public Pdiag_Double
{
	public:
	Parallel_Orbitals();
	~Parallel_Orbitals();

	// type : Sloc(1) Hloc(2) Hloc_fixed(3)
	bool in_this_processor(const int &iw1_all, const int &iw2_all);
	
	int* trace_loc_row;
	int* trace_loc_col;
	int out_hs; // mohan add 2010-09-02

	void set_trace(void);

		
};

#endif
