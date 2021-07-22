//==========================================================
// AUTHOR : mohan
// DATE : 2009-4-2
// Last Modify: 2009-08-28
//==========================================================
#ifndef NUMERICAL_BASIS_H
#define NUMERICAL_BASIS_H
#include "../src_pw/tools.h"
#include "bessel_basis.h"
//==========================================================
// CLASS :
// NAME :  Numerical_Basis 
//==========================================================
class Numerical_Basis
{
	public:
	Numerical_Basis();
	~Numerical_Basis();

	void start_from_file_k( const int &ik, ComplexMatrix &psi);
	void output_overlap( const ComplexMatrix *psi);

	private:

	static bool init_label;

	static Bessel_Basis bessel_basis;

	static IntArray *mu_index;
	static void init_mu_index(void);

	void numerical_atomic_wfc(const int &ik,const int &np,ComplexMatrix &psi);

	void Sq_overlap( 
		realArray &Sq_imag,
		realArray &Sq_real,
		const int &ik, 
		const int &np,
		const int derivative_order );

	void jlq3d_overlap(realArray &overlap_Q1, realArray &overlap_Q2,
		const int &ik_ibz, const int &ik, const int &np, const ComplexMatrix &psi,
		const int derivative_order );
		
	matrix psi_overlap(const ComplexMatrix *psi, const int derivative_order) const;
	
	void output_overlap_Sq(
		const string &name,
		ofstream &ofs, 
		const realArray *Sq_real, const realArray *Sq_imag);

	void output_overlap_Q(
		ofstream &ofs,
		const realArray &overlap_Q1,
		const realArray &overlap_Q2);

	void output_overlap_V(
		ofstream &ofs,
		const matrix &overlap_V) const;

};

#endif
