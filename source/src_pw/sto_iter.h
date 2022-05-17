#ifndef STO_ITER_H
#define STO_ITER_H
#include "sto_wf.h"
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

    void init(const int, int* nchip_in);
    
    void sum_stoband(Stochastic_WF& stowf);

    double calne(void);

    void itermu(int & iter);

    void sumpolyval_k(const int &ik, Stochastic_WF& stowf);

    void orthog(const int &ik, Stochastic_WF& stowf);

    void checkemm(const int &ik, int &iter, Stochastic_WF& stowf);

    Stochastic_Chebychev stoche;

    Stochastic_hchi stohchi;


	static double mu;
	static double mu0; // chemical potential; unit in Ry

    static double Emin;
	static double Emax; // unit in Ry

    bool change;
    
    //dos
    static double fwhm;
    static double targ_e;

    double targetne;

    double *spolyv;

	public:
    
    int * nchip;
    double th_ne;
    double KS_ne;

    public:
    static double root_fd(double e);
    static double fd(double e);
    static double nroot_fd(double e);
    static double nfd(double e);
    static double nxfd(double e);
    static double fdlnfd(double e);
    static double nfdlnfd(double e);

};

#endif// Eelectrons_Iter
