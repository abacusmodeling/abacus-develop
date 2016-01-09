#ifndef USE_HAMILT_MATRIX_H
#define USE_HAMILT_MATRIX_H

#include "../src_pw/tools.h"
#include "use_overlap_matrix.h"
#include "gint_gamma.h"
#include "gint_speed.h" //mohan add 2012-03-29
#include "gint_k.h"
#include "grid_integral.h"

class Use_Hamilt_Matrix
{
	public:
	Use_Hamilt_Matrix();
	~Use_Hamilt_Matrix();

    void set_ion(void);
		
	// used fro k-dependent Hamiltonian matrix.
	void calculate_Hk( const int &ik);
	
	// used for Gamma only Hamiltonian matrix.
	void calculate_Hgamma(void);

	// used for gamma only algorithms.
	Gint_Gamma GG;
	// fast grid integration, mohan add 2012-03-29
	Gint_Speed GS;

	// used for k-dependent grid integration.
	Gint_k GK;

	// use overlap matrix.
	Use_Overlap_Matrix UOM;

	// init S (overlap matrix) flag.
    bool init_s;

	private:

	// used for gamma only algorithms.
	void calculate_STNR_gamma(void);
	void calculate_STNR_gamma_B(void); //mohan add 2012-04-14
	void calculate_STNR_k(void);

};

#endif
