#ifndef SUBGRID_OPER_H
#define SUBGRID_OPER_H

#include <complex>
using namespace std;

class SubGrid_oper
{
	public:

	SubGrid_oper();
	~SubGrid_oper();

	//---------------------------
	// calculate the total wfc
	// occupations 
	//---------------------------

	//---------------------------
	// distribute the total wfc
	// into each proc of grid
	// group
	//---------------------------

	int* trace_lo_tot;
	double*** totwfc;
	complex<double>*** totwfc_B;//mohan add 2012-04-13
	bool allocate_totwfc;
	int lgd;

	private:


};
#endif
