#include "./math_chebyshev.h" 
#include "./constants.h"
#include "./blas_connector.h"
#include "./tool_quit.h"
#include "./global_function.h"
namespace ModuleBase
{
//we only have two examples: double and float.
template class Chebyshev<double>;
#ifdef __MIX_PRECISION
template class Chebyshev<float>;
#endif

// template<>
// FFTW<double>::FFTW(const int norder2_in)
// {
//     ccoef = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * norder2_in);
//     dcoef = (double *) fftw_malloc(sizeof(double) * norder2_in);
//     plancoef = fftw_plan_dft_r2c_1d(norder2_in, dcoef, ccoef, FFTW_ESTIMATE);
// }
// template<>
// FFTW<double>::~FFTW()
// {
//     fftw_destroy_plan(plancoef);
// 	fftw_free(ccoef);
// 	fftw_free(dcoef);
// }

// #ifdef __MIX_PRECISION
// template<>
// FFTW<float>::FFTW(const int norder2_in)
// {
//     ccoef = (fftwf_complex *) fftw_malloc(sizeof(fftwf_complex) * norder2_in);
//     dcoef = (float *) fftw_malloc(sizeof(float) * norder2_in);
//     plancoef = fftwf_plan_dft_r2c_1d(norder2_in, dcoef, ccoef, FFTW_ESTIMATE);
// }
// template<>
// FFTW<float>::~FFTW()
// {
//     fftwf_destroy_plan(plancoef);
// 	fftw_free(ccoef);
// 	fftw_free(dcoef);
// }
// #endif

//A number to control the number of grids in C_n integration
#define EXTEND 16

template<typename REAL>
Chebyshev<REAL>::Chebyshev(const int norder_in, const int ndmax_in) : fftw(2 * EXTEND * norder_in)
{
    this->norder = norder_in;
    norder2 = 2 * norder * EXTEND;
    if(this->norder < 1)
    {
        ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "The Chebyshev expansion order should be at least 1!");
    }
    polyvalue = new REAL [norder];
    coef_real = new REAL [norder];
    coef_complex = new std::complex<REAL> [norder];

    ndmin = ndmax = ndmax_in;

    getcoef = false;
    getpolyval = false;
}



template<typename REAL>
Chebyshev<REAL>::~Chebyshev()
{
	delete [] polyvalue;
	delete [] coef_real;
    delete [] coef_complex;
}

// REAL Chebyshev<REAL>::sumallterms(void)
// {
//     if(!getcoef||!getpolyval) 
// 	{
// 		ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef_real or polyval first!");
// 	}

//     REAL result = 0.0;
//     for(int ior = 0; ior < norder; ++ior)
//     {
//         result += coef_real[ior] * polyvalue[ior];
//     }
//     return result;   
// }
// template<typename REAL>
// void Chebyshev<REAL>:: calpolyval( REAL x )
// {
//     polyvalue[0] = 1;
//     polyvalue[1] = x;
//    for(int ior = 2; ior < this->norder; ++ior)
//    {
//        polyvalue[ior] = 2*x*polyvalue[ior-1]-polyvalue[ior-2];
//    }
//    getpolyval = true;
// }

template<typename REAL>
void Chebyshev<REAL>::calpolyval(
	void tfun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein, 
	const int m)
{
    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new std::complex<REAL> [ndmxt];
    arrayn = new std::complex<REAL> [ndmxt];
    arrayn_1 = new std::complex<REAL> [ndmxt];

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
        
        std::complex<REAL>* tem = arrayn_1;
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

template<typename REAL>
void Chebyshev<REAL>::calpolyvec(
	void tfun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein, std::complex<REAL> *polyvec,
	const int m)
{
    int disp;
    if(m==1) disp = this->ndmin;
    else
    {
        disp = this->ndmax * m;
    }
    
    std::complex<REAL> *arraynp1 = polyvec + 2 * disp;
	std::complex<REAL> *arrayn = polyvec + disp;
	std::complex<REAL> *arrayn_1 = polyvec;

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

template<typename REAL>
void Chebyshev<REAL>::calfinalvec_real(
	void tfun(std::complex<REAL> *in, std::complex<REAL> *out, const int n), 
	std::complex<REAL> *wavein, 
	std::complex<REAL> *waveout, 
	const int m)
{
    if(!getcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef first!");

    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new std::complex<REAL> [ndmxt];
    arrayn = new std::complex<REAL> [ndmxt];
    arrayn_1 = new std::complex<REAL> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    tfun(arrayn_1, arrayn,m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef_real[0] * arrayn_1[i] + coef_real[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef_real[ior] * arraynp1[i];
        }
        std::complex<REAL> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

template<typename REAL>
void Chebyshev<REAL>::calfinalvec_complex(
	void tfun(std::complex<REAL> *in, std::complex<REAL> *out, const int n), 
	std::complex<REAL> *wavein, 
	std::complex<REAL> *waveout, 
	const int m)
{
    if(!getcoef) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef_real first!");

    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;

    int ndmxt = ndmax * m;

    arraynp1 = new std::complex<REAL> [ndmxt];
    arrayn = new std::complex<REAL> [ndmxt];
    arrayn_1 = new std::complex<REAL> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    tfun(arrayn_1, arrayn,m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef_complex[0] * arrayn_1[i] + coef_complex[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs(arraynp1, arrayn, arrayn_1, tfun, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef_complex[ior] * arraynp1[i];
        }
        std::complex<REAL> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

template<typename REAL>
bool Chebyshev<REAL>::checkconverge(
	void tfun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein,
	REAL& tmax, 
	REAL &tmin, 
	REAL stept)
{
    bool converge = true;
    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;

    arraynp1 = new std::complex<REAL> [ndmin];
    arrayn = new std::complex<REAL> [ndmin];
    arrayn_1 = new std::complex<REAL> [ndmin];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmin);
    //LapackConnector::copy(ndim,wavein,1,arrayn_1,1); 
    if(tmin == tmax) 
	{
		tmax += stept;
	}

    tfun(arrayn_1, arrayn,1);
    REAL sum1,sum2;
    REAL t;

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
            arraynp1[i]=2.0*arraynp1[i]-arrayn_1[i];
        }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
  
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return converge;
}

template<typename REAL>
REAL Chebyshev<REAL>:: ddot_real(const int &m,
    const std::complex<REAL>* psi_L,
    const std::complex<REAL>* psi_R)
{
    REAL result = 0;
    if(ndmin == ndmax)
    {
        int dim2=2*ndmin*m;
        REAL *pL,*pR;
        pL=(REAL *)psi_L;
        pR=(REAL *)psi_R;
        result=BlasConnector::dot(dim2,pL,1,pR,1);
    }
    else
    {
        REAL *pL,*pR;
        pL=(REAL *)psi_L;
        pR=(REAL *)psi_R;
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
}