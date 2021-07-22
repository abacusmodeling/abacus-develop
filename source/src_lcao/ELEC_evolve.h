#ifndef ELEC_EVOLVE_H
#define ELEC_EVOLVE_H

#include "../src_pw/tools.h"
#include "LCAO_hamilt.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to evolve the electronic wave functions
// in TDDFT in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class ELEC_evolve
{

	friend class ELEC_scf;

	public:

	ELEC_evolve();
	~ELEC_evolve();

	// fuxiang add 2021-05-25

    static int tddft;
    static double td_dr2;
    static double td_dt;
    static double td_force_dt;
    static int td_val_elec_01;
    static int td_val_elec_02;
    static int td_val_elec_03;
    static int td_vext;
    static int td_vext_dire;
	static double td_timescale;
	static int td_vexttype;
	static int td_vextout;
	static int td_dipoleout;


	private:

	static void evolve_psi(const int &istep, LCAO_Hamilt &uhm, complex<double>*** wfc);
	void evolve_complex_matrix(const int &ik, complex<double> **c, complex<double> **c_init)const;
	void using_LAPACK_complex(const int &ik, complex<double> **c, complex<double> **c_init)const;
	void using_LAPACK_complex_2(const int &ik, complex<double>** c, complex<double>** c_init)const;

};

#endif
