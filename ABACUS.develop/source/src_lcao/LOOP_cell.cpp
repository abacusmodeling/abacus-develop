#include "LOOP_cell.h"
#include "LOOP_ions.h"

#include "src_lcao/dftu.h"   //Quxin add for DFT+U on 20201029

LOOP_cell::LOOP_cell(){}
LOOP_cell::~LOOP_cell(){}

void LOOP_cell::opt_cell()
{

	// Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		dftu.init(ucell, ParaO);
	}

	LOOP_ions ions;
	ions.opt_ions();
	
	return;
}

