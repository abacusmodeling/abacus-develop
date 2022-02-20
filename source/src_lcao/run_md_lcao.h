#ifndef RUN_MD_LCAO_H
#define RUN_MD_LCAO_H 

#include "../src_pw/charge_extra.h"
#include "module_orbital/ORB_control.h"

class Run_MD_LCAO
{

	public:

	Run_MD_LCAO();
	~Run_MD_LCAO();

	void opt_cell(ORB_control &orb_con);
	void opt_ions(void);
	void md_force_virial(const int &istep,
        const int& numIon, 
        double &potential, 
        ModuleBase::Vector3<double>* force, 
        ModuleBase::matrix& virial);

	private:

	// electron charge density extropolation method
	Charge_Extra CE;
	bool cellchange;
};

#endif
