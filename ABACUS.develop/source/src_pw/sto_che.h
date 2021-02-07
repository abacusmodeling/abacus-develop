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
    void calcoef(double fun(double));
    void recurs(double&tnp1, double &tn, double &tn_1, double& t); //tnp1: T_(n+1), tn: T_n, tn_1: T_(n-1)
    template<class T>
    void recurs(T *arraynp1, T* arrayn, T *arrayn_1, void fun(T *in,T *out), int& ndim);
    int norder;
    int norder2;  // 2 * norder
    double* coef;  //[norder] expansion coefficient of each order,
    complex<double> *ccoef;  //[norder2] temporary complex expansion coefficient of each order, only first norder coefficient is usefull.
    fftw_plan plancoef;
    bool initplan, initcoef;

	private:



};

template<class T>
void Stochastic_Chebychev:: recurs(T *arraynp1, T* arrayn, T *arrayn_1, void fun(T *in,T *out), int& ndim)
{
    fun(arrayn,arraynp1,ndim);
    for(int i = 0; i < ndim; ++i)
    {
        arraynp1[i]=2*arraynp1[i]-arrayn_1[i];
    }
}

#endif// Eelectrons_Chebychev
