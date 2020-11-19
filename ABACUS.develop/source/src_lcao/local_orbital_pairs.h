#ifndef LOCAL_ORBTIAL_PAIRS_H
#define LOCAL_ORBITAL_PARIS_H

#include "../src_pw/tools.h"

class Auxiliary_Basis
{
	public:

	
};

class Local_Orbital_Pairs
{
	public:
	Local_Orbital_Pairs();
	~Local_Orbital_Pairs();

	void init(void);

	private:
	int naux; // number of auxiliary basis.
	int lmax_aux; // DIY max L of auxiliary basis.
	
	double **aux_basis; // (naux, mesh) auxiliary basis.
	int mesh; // number of mesh points.
	double *r; // radial wave functions.
	double *rab; // rab
	

};

#endif
