#ifndef STOCHASTIC_WF_H
#define STOCHASTIC_WF_H

#include "tools.h"
#include "sto_che.h"
//qianrui 2021-2-4

//----------------------------------------------
// Generate stochastic wave functions 
//----------------------------------------------

class Stochastic_WF 
{
	public:

    Stochastic_WF();

    ~Stochastic_WF();


	void init();

	// ComplexMatrix may not be a best filetype to store the electronic wave functions
    ComplexMatrix* chi0;  	// origin stochastic wavefunctions in real space
	ComplexMatrix* chiortho;	// stochastic wavefunctions after in reciprocal space orthogonalized with KS wavefunctions
	Stochastic_Chebychev sto_che;
	int nchi; 				// Total number of stochatic obitals; unit in a_0^(3/2)
	int nchip; 				// The number of stochatic obitals in current process 
	int nche_sto;			// number of orders for Chebyshev expansion
	unsigned int seed_sto;  //random seed 

	int nbands_diag; // number of bands obtained from diagonalization

	int nbands_total; // number of bands in total, nbands_total=nchi+nbands_diag;

	protected:




	// check the properties of chi
    void check_chi(const ComplexMatrix *chi)const;

};
#endif 
