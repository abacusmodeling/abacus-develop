#ifndef INCLUDE_STO_ITER_H
#define INCLUDE_STO_ITER_H

#include "tools.h"
#include "sto_che.h"
#include "sto_hchi.h"

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

    void init();
    
    void calrho(double *);
    double calne();
    double caldnedmu();
    void itermu();
    void sumpolyval();
    Stochastic_Chebychev stoche;
    Stochastic_hchi stohchi;


	static double mu; // chemical potential
    static double Emin, Emax;
    double targetne;
    double * spolyv;

	private:
    
    double nchip;
    double th_ne;

    static double root_fd(double e);
    static double dfddmu(double e);
    static double fd(double e);
    static double nroot_fd(double e);
    static double ndfddmu(double e);
    static double nfd(double e);


};

#endif// Eelectrons_Iter
