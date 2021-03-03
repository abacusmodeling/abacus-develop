#ifndef INCLUDE_STO_CHEBYCHEV_H
#define INCLUDE_STO_CHEBYCHEV_H

#include "tools.h"
#include "../src_global/mymath.h"
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
    complex<double> calresult();
    void calresult(double &t, double &result);
    
    template<class T>
    void calresult(void fun(T *in, T *out), int& ndim, T *wavein, T *waveout);



    void calpolyval(void fun(complex<double> *in, complex<double> *out), int& ndim, complex<double> *wavein);

    int norder;
    int extend;
    int norder2;  // 2 * norder
    double* coef;  //[norder2] expansion coefficient of each order, only first norder coefficients are usefull
    double* dcoef;
    fftw_complex *ccoef;  //[norder2] temporary complex expansion coefficient of each order, only first norder coefficients are usefull.
    double *polyvalue; //
    fftw_plan plancoef;
    bool initplan, initcoef, getcoef, getpolyval;

	private:
    void recurs(double&tnp1, double &tn, double &tn_1, double& t); //tnp1: T_(n+1), tn: T_n, tn_1: T_(n-1)
    template<class T>
    void recurs(T *arraynp1, T* arrayn, T *arrayn_1, void fun(T *in,T *out), int& ndim);


};

template<class T>
void Stochastic_Chebychev:: recurs(T *arraynp1, T* arrayn, T *arrayn_1, void fun(T *in,T *out), int& ndim)
{
    fun(arrayn,arraynp1);
    for(int i = 0; i < ndim; ++i)
    {
        arraynp1[i]=2*arraynp1[i]-arrayn_1[i];
    }
}

template<class T>
void Stochastic_Chebychev:: calresult(void tfun(T *in, T *out), int &ndim, T *wavein, T *waveout)
{
    if(!getcoef) WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");

    T *arraynp1, *arrayn, *arrayn_1;
    arraynp1 = new T [ndim];
    arrayn = new T [ndim];
    arrayn_1 = new T [ndim];
    DCOPY(wavein, arrayn_1, ndim);
    tfun(arrayn_1, arrayn);
    
    //0- & 1-st order
    for(int i = 0; i < ndim; ++i)
    {
        waveout[i] = coef[0] * arrayn_1[i] + coef[1] * arrayn[i];
    }

    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, ndim);
        for(int i = 0; i < ndim; ++i)
        {
            waveout[i] += coef[ior] * arraynp1[i];
        }
        T * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}





#endif// Eelectrons_Chebychev
