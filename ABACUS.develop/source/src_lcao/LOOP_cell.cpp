#include "LOOP_cell.h"
#include "LOOP_ions.h"


RELAX_cell::RELAX_cell(){}

RELAX_cell::~RELAX_cell(){}


void RELAX_cell::opt_cell()
{
	LOOP_ions ions;
	ions.opt_ions();
	
	return;
}

