#ifndef GRID_BASE_BETA_H
#define GRID_BASE_BETA_H

#include "../src_pw/tools.h"
#include "../module_base/matrix3.h"
#include "../module_orbital/ORB_atomic_lm.h"

//AUTHOR : mohan
//DATE : 2008-09-16
// this class is inherited by Grid_Integral_Beta.h
// this class provides basic Grid operation and the 
// corresponding information.
class Grid_Base_Beta
{
	friend class Gint_Gamma;
	
public:

	void prepare( 
		const ModuleBase::Matrix3 &latvec_in, 
		const double &lat0_in);

protected:

	Grid_Base_Beta();
	~Grid_Base_Beta();
	
//==========================================================
// EXPLAIN : ModuleBase::Integral On 3D Real Space For Local Potential
// MEMBER FUNCTION :
//===========================================================
	double* rho1; // about charge
	double **density_kernel;
	double vfactor;
	ModuleBase::Matrix3 latvec;
	ModuleBase::Matrix3 latvec0;
	int* nnn;
	double lat0;
	enum cal_type{ cal_charge, cal_local, cal_vnlb } job;
};

#endif
