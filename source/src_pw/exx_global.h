#ifndef EXX_GLOBAL_H
#define EXX_GLOBAL_H

#include "xc_type.h"

struct Exx_Global
{
	enum class Hybrid_Type {No,HF,PBE0,HSE,Generate_Matrix};
	struct Exx_Info
	{
		Exx_Global::Hybrid_Type hybrid_type;

		double hybrid_alpha = 0.25;
		double hse_omega = 0.11;
		
		bool separate_loop = true;
		size_t hybrid_step = 1;
		
		void set_xcfunc( xcfunc &xcf ) const
		{
			switch( hybrid_type )
			{
				case Exx_Global::Hybrid_Type::HF:
					xcf.iexch_now=5;	xcf.igcx_now=0;
					xcf.icorr_now=0;	xcf.igcc_now=0;
					break;
				case Exx_Global::Hybrid_Type::PBE0:
					xcf.iexch_now=6;	xcf.igcx_now=8;
					break;
				case Exx_Global::Hybrid_Type::HSE:
					xcf.iexch_now=9;	xcf.igcx_now=12;
					break;
				default:
					throw std::invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			}			
		}			
	};
	Exx_Info info;
};

#endif
