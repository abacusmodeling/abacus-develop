#include "sto_che.h" 
#include "global.h"


Stochastic_Chebychev::Stochastic_Chebychev()
{
    initplan = false;
    norder = 0;
}

Stochastic_Chebychev::~Stochastic_Chebychev()
{
    fftw_destroy_plan(plancoef);
    delete [] coef;
    delete [] ccoef;
}
void Stochastic_Chebychev:: init()
{
    if(norder != 0)
    {
        norder2 = 2 * norder;
        coef = new double [norder2];
        ccoef = new complex<double> [norder2];
        initcoef = true;
    }
    else
    {
        WARNING_QUIT("Stochastic_Chebychev", "The Chebychev expansion order should be at least one!");
    }

}
    
void Stochastic_Chebychev:: calcoef(double* fun(double))
{
    if(!initcoef) WARNING_QUIT("Stochastic_Chebychev", "Please init coef first!");
    for(int i = 0; i < norder2; ++i)
    {
        coef = fun(cos((i+0.5)*PI/norder));
    }
     if(!initplan)
    {
        initplan = true;
        plancoef = fftw_plan_dft_1d(norder2, (fftw_complex *) coef, (fftw_complex *) ccoef, FFTW_FORWARD, FFTW_MEASURE);
    }
    fftw_execute(plancoef);
    for(int i = 0; i<norder2; ++i)
    {
        if(i == 0)
        {
            coef[i]=ccoef[i].real()/norder2;
        }
        else
        {
            coef[i]=ccoef[i].real()/norder;
        }
    }
}