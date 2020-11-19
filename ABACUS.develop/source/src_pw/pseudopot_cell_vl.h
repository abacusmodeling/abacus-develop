//==========================================================
// AUTHOR : Lixin He,mohan
// DATE : 2008-11-6
//==========================================================
#ifndef pseudopot_cell_vl_H
#define pseudopot_cell_vl_H

#include "tools.h"

class pseudopot_cell_vl
{
public:

	pseudopot_cell_vl();
	~pseudopot_cell_vl();

	matrix vloc;//(ntype,ngl),the local potential for each atom type(ntype,ngl)
	bool *numeric;//[ntype], =true

	void init_vloc(void);

private:

	// MODULE pseud
	// The variables describing pseudopotentials in analytical form
//	matrix cc;    // (2,npsx); the coefficients of the erf functions
//	matrix alpc;  // (2,npsx); the alpha of the erf functions
	double *zp;   // (npsx),the charge of the pseudopotential
//	realArray aps;   // (6,0:3,npsx);  the a_l coefficient
//	realArray alps;  // (3,0:3,npsx);  the b_l coefficient

//	double *a_nlcc;     // (npsx), nonlinear core correction coefficients:
//	double *b_nlcc;     // (npsx), rho_c(r) = (a_c + b_c*r^2) exp(-alpha_c*r^2)
//	double *alpha_nlcc; // (npsx)

//	int *nlc;          //(npsx), number of erf functions
//	int *nnl;          //(npsx), number of the gaussian functions
//	int *lmax;         //(npsx), maximum angular momentum of the pseudopot
//	int *lloc;         //(npsx), angular momentum of the part taken as local
	// END MODULE pseud

	void allocate(void);

	// generate vloc for a particular atom type.
	void vloc_of_g( const int &msh, const double *rab, const double *r, const double *vloc_at,
	               const double &zp, double *vloc) const;

	void print_vloc(void) const;

};

#endif // PSEUDOPOT_CELL_VL_H
