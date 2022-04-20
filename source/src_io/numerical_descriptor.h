//==========================================================
// AUTHOR : mohan
// DATE : 2021-01-04
//==========================================================
#ifndef NUMERICAL_DESCRIPTOR_H
#define NUMERICAL_DESCRIPTOR_H
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "bessel_basis.h"
//==========================================================
// CLASS :
// NAME :  Numerical_Descriptor 
//==========================================================
class Numerical_Descriptor
{
	public:
	Numerical_Descriptor();
	~Numerical_Descriptor();

	void output_descriptor( const ModuleBase::ComplexMatrix *psi, const int &lmax_in); // mohan added 2021-01-03

	private:

	bool init_label;

	int lmax; // lmax for descriptor
	int nmax; // nmax for descriptor
	int nlocal; // total number of atomic orbitals

	Bessel_Basis bessel_basis;

	ModuleBase::IntArray *mu_index;
	void init_mu_index(void);//mohan added 2021-01-03

	void jlq3d_overlap(ModuleBase::realArray &overlap_Q1, ModuleBase::realArray &overlap_Q2,
		const int &ik_ibz, const int &ik, const int &np, const ModuleBase::ComplexMatrix &psi);

	void generate_descriptor(ModuleBase::realArray &overlap_Q1, ModuleBase::realArray &overlap_Q2, 
		const int &it, const int &ia, double *d, const int &nd);

};

#endif
