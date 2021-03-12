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
    
    void sum_stoband();
    double calne();
    void itermu(int & iter);
    void sumpolyval();
    void orthog();
    void checkemm(int &iter);
    void test(); //only for test
    Stochastic_Chebychev stoche;
    Stochastic_hchi stohchi;


	static double mu, mu0; // chemical potential; unit in Ry
    static double Emin, Emax; // unit in Ry
    double targetne;
    double * spolyv;
    string stotype;

	private:
    
    double nchip;
    double th_ne;
    double KS_ne;

    static double root_fd(double e);
    static double fd(double e);
    static double nroot_fd(double e);
    static double nfd(double e);
    static double fdlnfd(double e);
    static double nfdlnfd(double e);


};

#endif// Eelectrons_Iter
