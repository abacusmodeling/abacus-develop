//=========================================================
//AUTHOR : liaochen, mohan
//DATE : 2008-03-28
//=========================================================
#ifndef GRID_INTEGRAL_H
#define GRID_INTEGRAL_H

#include "../src_pw/tools.h"
#include "grid_base.h"
//=========================================================
//CLASS  Grid_Integral
//Note : Integral On 3D Grids
//Feature : Matrix Elements Of Local Potential For 
//Numerical Orbitals
//=========================================================

class Grid_Integral : public Grid_Base
{
public:

	Grid_Integral(){};
	~Grid_Integral(){};

	void cal_rho(void);

	void cal_vlocal( const double* vlocal_in, 
		SparseMatrix &SM_in);

private:

	void deal_region_atoms(void);

	void deal_region_orbitals
	(
	 	const Vector3<double> &R1, //atomic position of the first atom
		const int &T1,
		const int &I1,
		const Vector3<double> &R2, //atomic position of the second atom
		const int &T2,
		const int &I2
	);

	void vlocal_in_small_box(void);
	void rho_in_small_box(void);


};

#endif
