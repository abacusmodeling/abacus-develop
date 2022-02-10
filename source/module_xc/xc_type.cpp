#include "xc_type.h"
#include "../module_base/global_function.h"
#include "../src_pw/global.h"
#include "exx_global.h"

xcfunc::xcfunc()
{
}

xcfunc::~xcfunc()
{
}


// Peize Lin add 2016-12-03
void xcfunc::hybrid_first()
{
#ifdef __LCAO
	// may do something
	ModuleBase::WARNING("functional","file "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__)+" may error ");
	if(Exx_Global::Hybrid_Type::HF==GlobalC::exx_global.info.hybrid_type)
	{
		iexch_now = 1;
		igcx_now = 3;
	}
	else if(Exx_Global::Hybrid_Type::PBE0==GlobalC::exx_global.info.hybrid_type)
	{
		iexch_now = 1;
		igcx_now = 3;
	}
	else if(Exx_Global::Hybrid_Type::HSE==GlobalC::exx_global.info.hybrid_type)
	{
		iexch_now = 1;
		igcx_now = 3;
	}
#endif
}
