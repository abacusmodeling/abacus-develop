#include "sto_che.h" 
#include "global.h"


Stochastic_Chebychev::Stochastic_Chebychev()
{
    initplan = false;
    initcoef = false;
    getcoef = false;
    getpolyval = false;
    norder = 10;
    ccoef = new complex<double> [1];
    coef = new double [1];
    polyvalue = new complex<double> [1];
}

Stochastic_Chebychev::~Stochastic_Chebychev()
{
    if(initplan)
    {
        fftw_destroy_plan(plancoef);
    }
    delete [] coef;
    delete [] ccoef;
    delete [] polyvalue;
}
void Stochastic_Chebychev:: init()
{
    norder = STO_WF.nche_sto;
    assert(norder > 10);
    if(norder != 0)
    {
        norder2 = 2 * norder;
        delete[] coef;
        delete[] ccoef;
        delete[] polyvalue;
        ccoef = new complex<double> [norder2];
        coef = new double [norder];
        polyvalue = new complex<double> [norder];
        ZEROS(polyvalue, norder);
        initcoef = true;
    }
    else
    {
        WARNING_QUIT("Stochastic_Chebychev", "The Chebychev expansion order should be at least one!");
    }

}
    
void Stochastic_Chebychev:: calcoef(double fun(double))
{
    if(!initcoef) WARNING_QUIT("Stochastic_Chebychev", "Please init coef first!");
    for(int i = 0; i < norder2; ++i)
    {
        ccoef[i]=complex<double>(fun(cos((i+0.5)*PI/norder)));
    }
     if(!initplan)
    {
        initplan = true;
        plancoef = fftw_plan_dft_1d(norder2, (fftw_complex *) ccoef, (fftw_complex *) ccoef, FFTW_FORWARD, FFTW_MEASURE);
    }
    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef[i]=ccoef[i].real()/norder2;
        }
        else
        {
            coef[i]/=ccoef[i].real()/norder;
        }
    }
    getcoef = true;
}

void Stochastic_Chebychev:: recurs(double&tnp1, double& tn, double& tn_1, double &t)
{   
    tnp1 = 2*t*tn-tn_1;
}

complex<double> Stochastic_Chebychev:: calresult()
{
    if(!getcoef||!getpolyval) WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef or polyval first!");
    complex<double> result = 0;
    for(int ior = 0; ior < norder; ++ior)
    {
        result += coef[ior] * polyvalue[ior];
    }
    return result;   
}

void Stochastic_Chebychev:: calresult(double &t, double& result)
{
    if(!getcoef) WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");
    double tnp1, tn, tn_1;
    tn_1 = 1;
    tn = tn_1 * t;
    //0- & 1-st order
    result = coef[0] * tn_1 + coef[1] * tn;
   
    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(tnp1, tn, tn_1, t);
        result += coef[ior] * tnp1;
        tn_1 = tn;
        tn = tnp1; 
    }
    return;
}

void Stochastic_Chebychev:: calpolyval(void tfun(complex<double> *in, complex<double> *out), int& ndim, complex<double> *wavein)
{
    if(!getcoef) WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");

    complex<double> *arraynp1, *arrayn, *arrayn_1;
    arraynp1 = new complex<double> [ndim];
    arrayn = new complex<double> [ndim];
    arrayn_1 = new complex<double> [ndim];

    for(int i = 0; i < ndim; ++i)
    {
        arrayn_1[i] = wavein[i]; 
    }
    tfun(arrayn_1, arrayn);

    //0- & 1-st order
    polyvalue[0] = ndim; // 0-th order : <wavein | wavein> = ndim
    for(int i = 0; i < ndim; ++i) // 1-st order : <wavein | tfun | wavein>
    {
        polyvalue[1] += conj(wavein[i]) * arrayn_1[i];
    }

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, ndim);
        for(int i = 0; i < ndim; ++i) // n-th order : <wavein | T_n(tfun) | wavein>
        {
            polyvalue[ior] += conj(wavein[i]) * arraynp1[i];
        }
        complex<double>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }

    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    getpolyval = true;
    return;
}
