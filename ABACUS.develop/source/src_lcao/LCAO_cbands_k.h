#ifndef LCAO_CBANDS_K_H
#define LCAO_CBANDS_K_H

#include "../src_pw/tools.h"
#include "use_hamilt_matrix.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to calculate the band structures
// (eigenvalues and eigen wave functions) of the 
// Kohn-Sham equation in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class LCAO_cbands_k
{

	friend class ELEC_scf;
	friend class ELEC_nscf;

	public:

	LCAO_cbands_k();
	~LCAO_cbands_k();


	private:

	static void cal_bands(const int &istep, Use_Hamilt_Matrix &uhm);


};

#endif
