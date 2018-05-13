//==========================================================
// AUTHOR : Lixin He,mohan
// DATE : 2008-11-08
//==========================================================
#ifndef PSEUDOPOT_CELL_VNL_H
#define PSEUDOPOT_CELL_VNL_H

#include "tools.h"
#include "pseudopot_cell_vl.h"

//==========================================================
// CLASS : 
// NAME : pseudopot_cell_vnl
// (Calculate the non-local pseudopotential in reciprocal
// space.
// Using plane wave as basis set.)
//==========================================================
class pseudopot_cell_vnl: public pseudopot_cell_vl
{
public:

//==========================================================
// MEMBER VARIABLES :
// NAME : dq(space between points in the pseudopotential tab)
// NAME : nkb(total number of beta functions, with struct.fact.)
// NAME : nqxq(size of interpolation table)
// NAME : nqx(number of interpolation points)
// NAME : nhm(max number of different beta functions per atom)
// NAME : lmaxkb(max angular momentum,(see pseudo_h))
// NAME : lmaxq(new added)
//==========================================================
	int nkb;
	int calculate_nqx(const double &ecutwfc,const double &dq);

	int nhm;
	int lmaxkb;
	int lmaxq;
//	int nbetam;		// max number of different projectors per atom ?
//  int nchim;		// max number of different wavefunctions per atom
//	int lllm;		// max number of

	matrix indv;		// indes linking  atomic beta's to beta's in the solid
	matrix nhtol;      	// correspondence n <-> angular momentum l
	matrix nhtolm;     	// correspondence n <-> combined lm index for (l,m)
	matrix nhtoj;		// new added

	realArray dvan;		//(:,:,:),  the D functions of the solid
	ComplexArray dvan_so;	//(:,:,:),  spin-orbit case,  added by zhengdy-soc
	realArray tab;		//(:,:,:), interpolation table for PPs
	realArray tab_at;	//(:,:,:), interpolation table for atomic wfc
	realArray deeq;		//(:,:,:,:), the integral of V_eff and Q_{nm}
	ComplexArray deeq_nc;	//(:,:,:,:), the spin-orbit case
	realArray becsum;	//(:,:,:,:), \sum_i  f(i) <psi(i)/beta_1><beta_m/psi(i)> //used in charge
//	realArray qq;		//(:,:,:), the q functions in the solid

	ComplexMatrix vkb;	// all beta functions in reciprocal space

	bool okvan;         // if .TRUE. at least one pseudo is Vanderbilt

	pseudopot_cell_vnl();
	~pseudopot_cell_vnl();

public:
	void init(const int ntype, const bool allocate_vkb=1);
	void init_vnl(void);
	void getvnl(const int &ik);

private:

	void print_vnl(ofstream &ofs);
};

#endif // PSEUDOPOT_CELL_VNL_H
