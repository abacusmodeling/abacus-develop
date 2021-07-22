#ifndef ELEC_CBANDS_K_H
#define ELEC_CBANDS_K_H

#include "../src_pw/tools.h"
#include "LCAO_hamilt.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to calculate the band structures
// (eigenvalues and eigen wave functions) of the 
// Kohn-Sham equation in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class ELEC_cbands_k
{

	friend class ELEC_scf;
	friend class ELEC_nscf;

	public:

	ELEC_cbands_k();
	~ELEC_cbands_k();


	private:

	static void cal_bands(const int &istep, LCAO_Hamilt &uhm);


};

#endif
