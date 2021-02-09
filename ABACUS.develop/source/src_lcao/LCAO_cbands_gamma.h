#ifndef LCAO_CBANDS_GAMMA_H
#define LCAO_CBANDS_GAMMA_H

#include "../src_pw/tools.h"
#include "use_hamilt_matrix.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to calculate the band structures
// (eigenvalues and eigen wave functions) of the 
// Kohn-Sham equation in terms of the gamma-only k-point
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class LCAO_cbands_gamma
{

	friend class Local_Orbital_Elec;

	public:

	LCAO_cbands_gamma();
	~LCAO_cbands_gamma();


	private:

	static void cal_bands(const int &istep, Use_Hamilt_Matrix &uhm);


};

#endif
