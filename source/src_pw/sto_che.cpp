#include "sto_che.h" 
#include "global.h"
#include "../module_base/constants.h"
#include "../module_base/lapack_connector.h"

Stochastic_Chebychev::Stochastic_Chebychev()
{
    initplan = false;
    initcoef = false;
    getcoef = false;
    getpolyval = false;
    extend = 16;
    norder = 5;
    polyvalue = new double [1];
    coef = new double [1];
    coefc = new complex<double>[1];
    ccoef = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *1);
    dcoef = (double *) fftw_malloc(sizeof(double) * 1);
    ndmin = ndmax = 0;
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
    delete [] coefc;
}

void Stochastic_Chebychev::init(int ndmax_in, int norder_in)
{
    this->ndmin = ndmax_in;
    this->ndmax = ndmax_in;
    this->norder = norder_in;

    if(norder<5)
    {
        ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "The Chebychev expansion order should be at least 5!");
    }

    assert(extend >= 1);
    fftw_free(ccoef);
    fftw_free(dcoef);
    delete [] polyvalue;
    delete [] coef;
    delete [] coefc;
    norder2 = 2 * norder * extend;
    ccoef = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * norder2);
    dcoef = (double *) fftw_malloc(sizeof(double) * norder2);
    plancoef = fftw_plan_dft_r2c_1d(norder2, dcoef, ccoef, FFTW_ESTIMATE);
    polyvalue = new double [norder];
    coef = new double [norder];
    coefc = new complex<double> [norder];
    initcoef = true;
    initplan = true;
    return;
}
    
void Stochastic_Chebychev::calcoef(double fun(double))
{
    if(!initcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please init coef first!");
    complex<double> *pcoef = (complex<double> *)ccoef;
    //three point = 2/3 M + 1/3 T;

    //(M)iddle point integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }

    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3;
        }
        else
        {
            coef[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3;
        }
    }

    //(T)rapezoid integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun(cos(i*ModuleBase::TWO_PI/norder2));
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
	return;
}

void Stochastic_Chebychev::calcoefc(double fun1(double),double fun2(double))
{
    if(!initcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please init coef first!");
    complex<double> *pcoef = (complex<double> *)ccoef;
    //three point = 2/3 M + 1/3 T;

    //(M)iddle point integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun1(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }

    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coefc[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coefc[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    //(T)rapezoid integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun1(cos(i*ModuleBase::TWO_PI/norder2));
    }
    
    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coefc[i].real( real(coefc[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coefc[i].real( real(coefc[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    //======================================================================

    //(M)iddle point integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun2(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }

    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coefc[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coefc[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    //(T)rapezoid integral method part
    for(int i = 0; i < norder2; ++i)
    {
        dcoef[i]=fun2(cos(i*ModuleBase::TWO_PI/norder2));
    }
    
    fftw_execute(plancoef);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coefc[i].imag( imag(coefc[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coefc[i].imag( imag(coefc[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 
    //==================================================================

    getcoef = true;
	return;
}


double Stochastic_Chebychev::sumallterms(void)
{
    if(!getcoef||!getpolyval) 
	{
		ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef or polyval first!");
	}

    double result = 0;
    for(int ior = 0; ior < norder; ++ior)
    {
        result += coef[ior] * polyvalue[ior];
    }
    return result;   
}

 void Stochastic_Chebychev:: calpolyval( double x )
 {
     polyvalue[0] = 1;
     polyvalue[1] = x;
    for(int ior = 2; ior < norder; ++ior)
    {
        polyvalue[ior] = 2*x*polyvalue[ior-1]-polyvalue[ior-2];
    }
    getpolyval = true;

 }

void Stochastic_Chebychev::calpolyval(
	void tfun(complex<double> *in, complex<double> *out, const int), 
	complex<double> *wavein, 
	const int m)
{
    complex<double> *arraynp1;
	complex<double> *arrayn;
	complex<double> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new complex<double> [ndmxt];
    arrayn = new complex<double> [ndmxt];
    arrayn_1 = new complex<double> [ndmxt];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);


    tfun(arrayn_1, arrayn,m);

    polyvalue[0] = this->ddot_real(m,wavein,wavein);
    polyvalue[1] = this->ddot_real(m,wavein,arrayn);
    
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
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        polyvalue[ior] = this->ddot_real(m,wavein,arraynp1);
        
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

void Stochastic_Chebychev::calpolyvec(
	void tfun(complex<double> *in, complex<double> *out, const int), 
	complex<double> *wavein, complex<double> *polyvec,
	const int m)
{
    int disp;
    if(m==1) disp = this->ndmin;
    else
    {
        disp = this->ndmax * m;
    }
    
    complex<double> *arraynp1 = polyvec + 2 * disp;
	complex<double> *arrayn = polyvec + disp;
	complex<double> *arrayn_1 = polyvec;

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, disp);

    tfun(arrayn_1, arrayn,m);

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        arrayn += disp;
        arrayn_1 += disp;
        arraynp1 += disp;
    }
    return;
}

void Stochastic_Chebychev::calfinalvec(
	void tfun(complex<double> *in, complex<double> *out, const int n), 
	complex<double> *wavein, 
	complex<double> *waveout, 
	const int m)
{
    if(!getcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");

    complex<double> *arraynp1;
	complex<double> *arrayn;
	complex<double> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new complex<double> [ndmxt];
    arrayn = new complex<double> [ndmxt];
    arrayn_1 = new complex<double> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    tfun(arrayn_1, arrayn,m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef[0] * arrayn_1[i] + coef[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef[ior] * arraynp1[i];
        }
        complex<double> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

void Stochastic_Chebychev::calfinalvec2(
	void tfun(complex<double> *in, complex<double> *out, const int n), 
	complex<double> *wavein, 
	complex<double> *waveout, 
	const int m)
{
    if(!getcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");

    complex<double> *arraynp1;
	complex<double> *arrayn;
	complex<double> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new complex<double> [ndmxt];
    arrayn = new complex<double> [ndmxt];
    arrayn_1 = new complex<double> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    tfun(arrayn_1, arrayn,m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coefc[0] * arrayn_1[i] + coefc[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coefc[ior] * arraynp1[i];
        }
        complex<double> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

void Stochastic_Chebychev::calfinalvec_complex(
	void tfun(complex<double> *in, complex<double> *out, const int n), 
	complex<double> *wavein, 
	complex<double> *waveout, 
	const int m)
{
    if(!getcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");

    complex<double> *arraynp1;
	complex<double> *arrayn;
	complex<double> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new complex<double> [ndmxt];
    arrayn = new complex<double> [ndmxt];
    arrayn_1 = new complex<double> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    tfun(arrayn_1, arrayn,m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coefc[0] * arrayn_1[i] + coefc[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coefc[ior] * arraynp1[i];
        }
        complex<double> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}


bool Stochastic_Chebychev::checkconverge(
	void tfun(complex<double> *in, complex<double> *out, const int), 
	complex<double> *wavein,
	double& tmax, 
	double &tmin, 
	double stept)
{
    bool converge = true;
    complex<double> *arraynp1;
	complex<double> *arrayn;
	complex<double> *arrayn_1;

    arraynp1 = new complex<double> [ndmin];
    arrayn = new complex<double> [ndmin];
    arrayn_1 = new complex<double> [ndmin];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmin);
    //LapackConnector::copy(ndim,wavein,1,arrayn_1,1); 
    if(tmin == tmax) 
	{
		tmax += stept;
	}

    tfun(arrayn_1, arrayn,1);
    double sum1,sum2;
    double t;

    sum1=ModuleBase::GlobalFunc::ddot_real(ndmin,arrayn_1,arrayn_1);
    sum2=ModuleBase::GlobalFunc::ddot_real(ndmin,arrayn_1,arrayn);
    t = sum2 / sum1 * (tmax - tmin) / 2 + (tmax + tmin) / 2;
    if(t < tmin || tmin == 0)
    {
        converge = false;
        tmin = t-stept;
    }
    if(t > tmax)
    {
        converge = false;
        tmax = t+stept;
    }

    for(int ior = 2; ior < norder; ++ior)
    {
        tfun(arrayn,arraynp1,1);
        sum1=ModuleBase::GlobalFunc::ddot_real(ndmin,arrayn,arrayn);
        sum2=ModuleBase::GlobalFunc::ddot_real(ndmin,arrayn,arraynp1);
        t = sum2/sum1 * (tmax - tmin) / 2 + (tmax + tmin) / 2;
        if(t < tmin)
        {
            converge = false;
            tmin = t - stept;
        }
        else if(t > tmax)
        {
            converge = false;
            tmax = t + stept;
        }
        for(int i = 0; i < ndmin; ++i)
        {
            arraynp1[i]=2*arraynp1[i]-arrayn_1[i];
        }
        complex<double>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
  
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return converge;
}
double Stochastic_Chebychev:: ddot_real(const int &m,
    const complex<double>* psi_L,
    const complex<double>* psi_R)
{
    double result = 0;
    if(ndmin == ndmax)
    {
        int dim2=2*ndmin*m;
        double *pL,*pR;
        pL=(double *)psi_L;
        pR=(double *)psi_R;
        result=BlasConnector::dot(dim2,pL,1,pR,1);
    }
    else
    {
        double *pL,*pR;
        pL=(double *)psi_L;
        pR=(double *)psi_R;
        for(int i = 0 ; i < m ; ++i)
        {
            int dim2=2*ndmin;
            result +=  BlasConnector::dot(dim2,pL,1,pR,1);
            pL += 2*ndmax;
            pR += 2*ndmax;
        }
    }
    return result;
}