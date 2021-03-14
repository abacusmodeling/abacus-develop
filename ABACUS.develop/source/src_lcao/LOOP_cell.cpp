#include "LOOP_cell.h"
#include "LOOP_ions.h"

#include "src_lcao/dftu.h"   //Quxin add for DFT+U on 20201029

LOOP_cell::LOOP_cell(){}
LOOP_cell::~LOOP_cell(){}

void LOOP_cell::opt_cell(void)
{
	TITLE("LOOP_cell","opt_cell");

	// Peize Lin 2016-12-03
	if (CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
	{
		switch(exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_lcao.init();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		}
	}	

	// Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		dftu.init(ucell, ParaO);
	}

	LOOP_ions ions;
	ions.opt_ions();
	
	return;
}

