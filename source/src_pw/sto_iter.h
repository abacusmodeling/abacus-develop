#ifndef STO_ITER_H
#define STO_ITER_H
#include "sto_wf.h"
#include "../module_base/math_chebyshev.h"
#include "sto_hchi.h"
#include "sto_func.h"

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

    ModuleBase::Chebyshev<double>* p_che;

    Stochastic_hchi stohchi;
    Sto_Func<double> stofunc;

	double mu0; // chemical potential; unit in Ry
    bool change;
    double targetne;
    double *spolyv;

	public:
    
    int * nchip;
    double th_ne;
    double KS_ne;

};

#endif// Eelectrons_Iter
