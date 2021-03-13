#include "LOOP_cell.h"
#include "LOOP_ions.h"


LOOP_cell::LOOP_cell(){}
LOOP_cell::~LOOP_cell(){}

void LOOP_cell::opt_cell()
{
	LOOP_ions ions;
	ions.opt_ions();
	
	return;
}

