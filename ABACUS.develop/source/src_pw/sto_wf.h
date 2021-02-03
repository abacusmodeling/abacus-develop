#ifndef STOCHASTIC_WF_H
#define STOCHASTIC_WF_H

#include "tools.h"

//----------------------------------------------
// Generate stochastic wave functions 
//----------------------------------------------

class Stochastic_WF 
{
	public:

    Stochastic_WF();

    ~Stochastic_WF();


	void allocate_chi();


	// ComplexMatrix may not be a best filetype to store the electronic wave functions
    ComplexMatrix *chi;  // stochastic wavefunctions in the PW basis

	int nchi;

	int nbands_diag; // number of bands obtained from diagonalization

	int nbands_total; // number of bands in total, nbands_total=nchi+nbands_diag;

	protected:


    // Calculate random wave functions as stochastic wave functions
    void random(ComplexMatrix &chi);

	// check the properties of chi
    void check_chi(const ComplexMatrix *chi)const;

};
#endif 
