#ifndef LCAO_EVOLVE_H
#define LCAO_EVOLVE_H

#include "../src_pw/tools.h"
#include "use_hamilt_matrix.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to evolve the electronic wave functions
// in TDDFT in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class LCAO_evolve
{

	friend class Local_Orbital_Elec;

	public:

	LCAO_evolve();
	~LCAO_evolve();

	private:

	static void evolve_psi(const int &istep, Use_Hamilt_Matrix &uhm, complex<double>*** wfc);


};

#endif
