#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "../src_pw/tools.h"
#include "FORCE_k.h"
//#include "./force_lcao_gamma.h"
#include "../src_pw/stress_func.h"
#include "../input_conv.h"
#include "../src_pw/forces.h"

class Force_Stress_LCAO  
{
	// mohan add 2021-02-09
	friend class md;
	friend class Run_MD_LCAO;
	friend void Input_Conv::Convert();
	friend class Update_input;
	friend class LOOP_ions;
	friend class MD_func;

	public :
	
	Force_Stress_LCAO ();
	~Force_Stress_LCAO ();

	private:

	Force_LCAO_k flk;
//	Force_LCAO_gamma flg;
	Stress_Func sc_pw;
	Forces f_pw;

	void allocate (void);
	void destroy (void);
	
	void print_force(const std::string &name, ModuleBase::matrix& f, const bool screen, bool ry)const;
	void printforce_total (const bool ry, const bool istestf, ModuleBase::matrix& fcs);
	
	void getForceStress(
		const bool isforce, 
		const bool isstress, 
		const bool istestf, 
		const bool istests, 
		ModuleBase::matrix &fcs, 
		ModuleBase::matrix &scs);

	void forceSymmetry(ModuleBase::matrix &fcs);

	void calForcePwPart(
		ModuleBase::matrix &fvl_dvl, 
		ModuleBase::matrix &fewalds, 
		ModuleBase::matrix &fcc, 
		ModuleBase::matrix &fscc);

	void calForceStressIntegralPart(
		const bool isGammaOnly,
		const bool isforce,
		const bool isstress,
		ModuleBase::matrix &foverlap,
		ModuleBase::matrix &ftvnl_dphi,
		ModuleBase::matrix &fvnl_dbeta,	
		ModuleBase::matrix &fvl_dphi,
		ModuleBase::matrix &soverlap,
		ModuleBase::matrix &stvnl_dphi,
		ModuleBase::matrix &svnl_dbeta,
		ModuleBase::matrix &svl_dphi);

	void calStressPwPart(
		ModuleBase::matrix &sigmadvl,
		ModuleBase::matrix &sigmahar,
		ModuleBase::matrix &sigmaewa,
		ModuleBase::matrix &sigmacc,
		ModuleBase::matrix &sigmaxc);
	
	static double force_invalid_threshold_ev;
	static double output_acc; // control the accuracy
};
#endif
