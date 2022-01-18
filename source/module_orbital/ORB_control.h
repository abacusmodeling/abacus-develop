#ifndef ORB_CONTROL_H 
#define ORB_CONTROL_H 

#include "ORB_gen_tables.h"
#include "ORB_read.h"

class ORB_control 
{
	public:

    ORB_control();
    ~ORB_control();

	//first step: read orbital file
	void read_orb_first(
		std::ofstream &ofs_in, 
		LCAO_Orbitals &orb,
		const int &ntype, // mohan add 2021-04-26
		const int &lmax, // mohan add 2021-04-26 
		const double &lcao_ecut_in, // mohan add 2021-04-16
		const double &lcao_dk_in, // mohan add 2021-04-16
		const double &lcao_dr_in, // mohan add 2021-04-16
		const double &lcao_rmax_in, // mohan add 2021-04-16
		const bool &deepks_setorb,
		const int &out_r_matrix,
		const bool &force_flag, // mohan add 2021-05-07
		const int &my_rank // mohan add 2021-04-26
	);
    /// Generate the S(overlap),T,NL matrix.
    void set_orb_tables(
		std::ofstream &ofs_in,
		ORB_gen_tables &OGT, 
		LCAO_Orbitals &orb,
		const double &lat0,
		const bool &deepks_setorb,
		const int &Lmax_exx,
		const int &nprojmax, 
		const int* nproj,
		const Numerical_Nonlocal* beta_); 

    void clear_after_ions(
		ORB_gen_tables &OGT, 
		LCAO_Orbitals &orb,
		const bool &deepks_setorb,
		const int* nproj_);

};
#endif
