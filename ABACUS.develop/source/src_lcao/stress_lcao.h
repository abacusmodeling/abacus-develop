#ifndef STRESS_LCAO_H
#define STRESS_LCAO_H
#include "./../src_pw/stress_func.h"

class Stress_LCAO:public Stress_Func
{
	public :
	
	Stress_LCAO (){};
	~Stress_LCAO (){};

	//calculate the stress in LCAO basis
	void start_stress(
		const double overlap[][3], 
		const double tvnl_dphi[][3],
		const double vnl_dbeta[][3],
		const double vl_dphi[][3], 
		const matrix& stress_vdw,
		matrix& scs);
};
#endif
