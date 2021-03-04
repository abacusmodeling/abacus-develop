#include "sto_che.h" 
#include "global.h"



Stochastic_Chebychev::Stochastic_Chebychev()
{
    initplan = false;
    initcoef = false;
    getcoef = false;
    getpolyval = false;
    extend = 2048;
    norder = 5;
}

Stochastic_Chebychev::~Stochastic_Chebychev()
{
    if(initplan)
    {
        fftw_destroy_plan(plancoef);
    }
    fftw_free(ccoef);
    fftw_free(dcoef);
    delete [] polyvalue;
    delete [] coef;
}
void Stochastic_Chebychev:: init()
{
    norder = STO_WF.nche_sto;
    assert(norder > 5);
    assert(extend >= 1);
    if(norder != 0)
    {
        norder2 = 2 * norder * extend;
        ccoef = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * norder2);//new complex<double> [norder2];
        dcoef = (double *) fftw_malloc(sizeof(double) * norder2);
        plancoef = fftw_plan_dft_r2c_1d(norder2, dcoef, ccoef, FFTW_ESTIMATE);
        polyvalue = new double [norder];
        coef = new double [norder];
        initcoef = true;
        initplan = true;
    }
    else
    {
        WARNING_QUIT("Stochastic_Chebychev", "The Chebychev expansion order should be at least one!");
    }

}
    
void Stochastic_Chebychev:: calcoef(double fun(double))
{
    if(!initcoef) WARNING_QUIT("Stochastic_Chebychev", "Please init coef first!");
    complex<double> *pcoef = (complex<double> *)ccoef;
    //three point = 2/3 M + 1/2 T;

    //(M)iddle point integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun(cos((i+0.5)*TWO_PI/norder2));
    }
    fftw_execute(plancoef);
    complex<double> ui(0,1);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef[i] = real(exp(ui*i*PI/norder2) * pcoef[i]) / norder2 * 2 / 3;
        }
        else
        {
            coef[i] = real(exp(ui*i*PI/norder2) * pcoef[i]) / norder2 * 4 / 3;
        }
    }

    //(T)rapezoid integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun(cos(i*TWO_PI/norder2));
    }
    
    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef[i] += real(pcoef[i]) / norder2 * 1 / 3;
        }
        else
        {
            coef[i] += real(pcoef[i]) / norder2 * 2 / 3;
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

    complex<double> *arraynp1, *arrayn, *arrayn_1;
    arraynp1 = new complex<double> [ndim];
    arrayn = new complex<double> [ndim];
    arrayn_1 = new complex<double> [ndim];

    //DCOPY(wavein, arrayn_1, ndim);

    LapackConnector::copy(ndim,wavein,1,arrayn_1,1);
    tfun(arrayn_1, arrayn);

    int inc =1 ;
    complex<double> overlap;
    zdotc_(&overlap,&ndim,wavein,&inc,wavein,&inc);
    polyvalue[0] = overlap.real();
    zdotc_(&overlap,&ndim,wavein,&inc,arrayn,&inc);
    polyvalue[1] = overlap.real();  
    //ZEROS(polyvalue,norder);
    ////0- & 1-st order
    //for(int i = 0; i < ndim; ++i) 
    //{
    //    polyvalue[0] += real(conj(wavein[i]) * wavein[i]);// 0-th order : <wavein | wavein> = ndim
    //    polyvalue[1] += real(conj(wavein[i]) * arrayn[i]);// 1-st order : <wavein | tfun | wavein>
    //}

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, ndim);
        zdotc_(&overlap,&ndim,wavein,&inc,arraynp1,&inc);
        polyvalue[ior] = overlap.real();
        //for(int i = 0; i < ndim; ++i) // n-th order : <wavein | T_n(tfun) | wavein>
        //{
        //    polyvalue[ior] += real(conj(wavein[i]) * arraynp1[i]);
        //}
        
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
