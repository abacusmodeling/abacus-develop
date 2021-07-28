//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_H
#define DFTU_H

#include <string>

#include "../src_pw/charge_broyden.h"
#include "dftu_relax.h"

#include "../module_cell/unitcell_pseudo.h"
#include "../src_parallel/parallel_orbitals.h"

using namespace std;

//==========================================================
// CLASS :
// NAME : DTFU (DFT+U)
//==========================================================
class DFTU : public DFTU_RELAX
{

public:

    DFTU();                      // constructor 
    ~DFTU();                     // deconstructor

	// initialize the input terms of  U, J, double_counting etc
    void init(UnitCell_pseudo &cell, // unitcell class
		Parallel_Orbitals &po // parallel orbitals parameters
	);
    
    void folding_overlap_matrix(const int ik, complex<double>* Sk);
    
    //calculate the local occupation number matrix
    void cal_occup_m_k(const int iter);
    void cal_occup_m_gamma(const int iter);

    void write_occup_m(const string &fn);
    void read_occup_m(const string &fn);
    void local_occup_bcast();
    
    //calculate the energy correction: en_cor
    //called at energy::calculate_etot(void)
    void cal_energy_correction( const int istep);

    //calculate the effective potential
    void cal_eff_pot_mat_complex(const int ik, const int istep, complex<double>* eff_pot);
    void cal_eff_pot_mat_real(const int ik, const int istep, double* eff_pot);
    void cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR);
    void cal_eff_pot_mat_R_complex_double(const int ispin, complex<double>* SR, complex<double>* HR);

    void print(const int T, const int iat, const int L, const int N, const int iter);

    void output();

    void out_orbitals();

    double EU;
    int iter_dftu;
   
};

extern DFTU dftu;

#endif
