#ifndef LCAO_MATRIX_H
#define LCAO_MATRIX_H

#include "../src_pw/tools.h"

class LCAO_Matrix 
{
	public:
	LCAO_Matrix();
	~LCAO_Matrix();

	//------------------------------
	// H, S, Hfixed 
	// used in gamma only algorithm.
	// thse matrix are used to
	// diagonalize.
	//------------------------------
	double* Hloc;
	double* Sloc;
	double* Hloc_fixed;
	
	//------------------------------
	// 1. Hamiltonian(vl),
	// 2. overlap matrix Sloc2
	// 3. fixed (vna+T+Vnl) matrix.
	// used in kpoint algorithm.
	// these matrix are used to
	// diagonalize.
	//------------------------------
	complex<double> *Hloc2;
	complex<double> *Sloc2;
	complex<double> *Hloc_fixed2;
	
	//------------------------------
	// Store H(mu,nu')
	// nu' : nu in near unitcell R.
	// used in kpoint algorithm.
	// these matrixed are used
	// for 'folding_matrix' in lcao_nnr,
	// HlocR -> Hloc2,
	// SlocR -> Sloc2, 
	//------------------------------
	double* HlocR;
	double* SlocR;
	double* Hloc_fixedR;


	//========================================
	// FORCE
	//========================================

	//-----------------------------------------
	// force in LCAO
	// used in gamma only algorithm.
	//-----------------------------------------
	double* DSloc_x;
	double* DSloc_y;
	double* DSloc_z;

	//-----------------------------------------
	// force in LCAO
	// used in k-points algorithm.
	//-----------------------------------------
	double* DSloc_Rx;
	double* DSloc_Ry;
	double* DSloc_Rz;

	//-----------------------------------------
	// dT + part of dVNL
	// used in gamma only algorithm.
	//-----------------------------------------
	double* DHloc_fixed_x; 
	double* DHloc_fixed_y;
	double* DHloc_fixed_z;

	//-----------------------------------------
	// dT + part of dVNL
	// used in kpoint algorithm.
	//-----------------------------------------
	double* DHloc_fixedR_x;
	double* DHloc_fixedR_y;
	double* DHloc_fixedR_z;


	void divide_HS_in_frag(void);
	void set_HSgamma(const int &iw1_all, const int &iw2_all, const double &v, const char &dtype);
	void set_HSk(const int &iw1_all, const int &iw2_all, const complex<double> &v, const char &dtype);
	void set_force (const int& iw1_all, const int& iw2_all, const double& vx, const double& vy, 
		const double& vz, const char &dtype);

	void zeros_HSgamma(const char &mtype);
	void zeros_HSk(const char &mtype);
	void zeros_HSR(const char &mtype, const int &nnr);

	void print_HSgamma(const char &mtype);
	void print_HSk(const char &mtype, const char &vtype = 'C', const double &accuracy = 1.0e-5);
	void update_Hloc(void);
	void update_Hloc2(void);

	void allocate_HS_R(const int &nnr);
	void allocate_HS_gamma(const int &nloc);
	void allocate_HS_k(const int &nloc);

	void output_HSk(const char &mtype, string &fn);

};

#endif
