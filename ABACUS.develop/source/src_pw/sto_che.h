#ifndef INCLUDE_STO_CHEBYCHEV_H
#define INCLUDE_STO_CHEBYCHEV_H

#include "tools.h"
#include "mymath.h"
//----------------------------------------------
// Chebychev Filtering
//----------------------------------------------

class Stochastic_Chebychev
{

	public:

    // constructor and deconstructor
    Stochastic_Chebychev();
    ~Stochastic_Chebychev();
    void init();
    void calcoef(double* fun(double));
    
    int norder;
    int norder2;  // 2 * norder
    complex<double> *ccoef; //temperary complex coefficient
    double *coef;  //expansion coefficient of each order, only first norder coefficient will be used.
    fftw_plan plancoef;
    bool initplan, initcoef;

	private:



};

#endif// Eelectrons_Chebychev
