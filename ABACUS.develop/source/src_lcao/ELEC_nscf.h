#ifndef ELEC_NSCF_H
#define ELEC_NSCF_H

#include "../src_pw/tools.h"
#include "use_hamilt_matrix.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to run non-self-consistent calculations 
// to solve the Kohn-Sham equation by reading in electron
// charge density
//-----------------------------------------------------------

class ELEC_nscf
{

	friend class Local_Orbital_Ions;
	friend class Local_Orbital_Elec;

	public:

	ELEC_nscf();
	~ELEC_nscf();

	private:

	static void nscf(Use_Hamilt_Matrix &uhm); 


};

#endif
