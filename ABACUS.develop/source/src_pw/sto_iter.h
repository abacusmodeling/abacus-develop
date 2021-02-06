#ifndef INCLUDE_STO_ITER_H
#define INCLUDE_STO_ITER_H

#include "tools.h"

//----------------------------------------------
// Solve for the new electron density and iterate 
// until SCF loop converges (mu, epsilon, rho)
// mu: chemical potential
// epsilon: eigenvalues
// rho: charge density
//----------------------------------------------

class Stochastic_Iter
{

	public:

    // constructor and deconstructor
    Stochastic_Iter();

    ~Stochastic_Iter();

	double mu; // chemical potential

	private:


};

#endif// Eelectrons_Iter
