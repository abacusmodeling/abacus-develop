#ifndef STO_CHEBYCHEV_H
#define STO_CHEBYCHEV_H
#include <complex>
#include "fftw3.h"

namespace ModuleBase
{
//template class for fftw
template<typename T>
class FFTW;

/**
 * @brief A class to treat the Chebyshev expansion.
 *
 * @author qianrui on 2022-05-18
 * @details
 * Math:
 * 
 * USAGE：
 *   
 *
 */
template<typename REAL>
class Chebyshev
{

public:

    // constructor and deconstructor
    Chebyshev(const int norder, const int ndmax_in);
    ~Chebyshev();

public:
    
	//Calculate coefficients C_n[f], where f is a function of real number
	template<class T>
    void calcoef_real(T *ptr, REAL (T::*fun)(REAL));
	//Calculate coefficients C_n[g], where g is a function of complex number
	template<class T>
    void calcoef_complex(T *ptr, std::complex<REAL> (T::*fun)(std::complex<REAL>));
	//Calculate coefficients C_n[g], where g is a general complex function g(x)=(g1(x), g2(x)) e.g. exp(ix)=(cos(x), sin(x))
	template<class T>
	void calcoef_pair(T *ptr, REAL (T::*fun1)(REAL), REAL (T::*fun2)(REAL));
    
    void calfinalvec_real(
		void fun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, 
		std::complex<REAL> *waveout, 
		const int m = 1);

	void calfinalvec_complex(
		void fun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, 
		std::complex<REAL> *waveout, 
		const int m = 1);

    bool checkconverge(
		void tfun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
 		std::complex<REAL> *wavein,
		REAL& tmax, 
		REAL& tmin, 
		REAL stept);

    void calpolyval(
		void fun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, 
		const int m =1);
	void calpolyvec(
		void fun(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, std::complex<REAL> *polyvec,
		const int m =1);

	// void calpolyval(REAL x);
	// REAL sumallterms();
	

    int norder;
    int norder2;  // 2 * norder * EXTEND

	// expansion coefficient of each order
	// only first norder coefficients are usefull
    REAL* coef_real; 
	std::complex<REAL>* coef_complex; 
	FFTW<REAL> fftw;

    REAL *polyvalue;

	bool getcoef;
	bool getpolyval;
	int ndmin; //dim of vector
	int ndmax; //the distance between two closed vectors

private:

    REAL ddot_real(const int &m,
    const std::complex<REAL>* psi_L,
    const std::complex<REAL>* psi_R);

    template<class T>
    void recurs(
		T *arraynp1, 
		T* arrayn, 
		T *arrayn_1, 
		void fun(T *in,T *out, const int),
		const int m);
};

template<typename REAL>
template<class T>
void Chebyshev<REAL>::recurs(
		T *arraynp1, 
		T* arrayn, 
		T *arrayn_1, 
		void fun(T *in,T *out, const int),
		const int m)
{
    fun(arrayn,arraynp1,m);
	for(int ib = 0 ; ib < m ; ++ib)
	{
    	for(int i = 0; i < ndmin; ++i)
    	{
        	arraynp1[i+ib*ndmax]=2.0*arraynp1[i+ib*ndmax]-arrayn_1[i+ib*ndmax];
    	}
	}
}

template<>
class FFTW<double>
{
public:
	FFTW(const int norder2_in)
	{
	    ccoef = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * norder2_in);
	    dcoef = (double *) fftw_malloc(sizeof(double) * norder2_in);
	    coef_plan = fftw_plan_dft_r2c_1d(norder2_in, dcoef, ccoef, FFTW_ESTIMATE);
	}
	~FFTW()
	{
		fftw_destroy_plan(coef_plan);
		fftw_free(ccoef);
		fftw_free(dcoef);
	}
    double* dcoef; //[norder2]
	fftw_complex *ccoef;
	fftw_plan coef_plan;
};

#ifdef __MIX_PRECISION
template<>
class FFTW<float>
{
public:
	FFTW(const int norder2_in)
	{
	    ccoef = (fftwf_complex *) fftw_malloc(sizeof(fftwf_complex) * norder2_in);
	    dcoef = (float *) fftw_malloc(sizeof(float) * norder2_in);
	    coef_plan = fftwf_plan_dft_r2c_1d(norder2_in, dcoef, ccoef, FFTW_ESTIMATE);
	}
	~FFTW()
	{
		fftwf_destroy_plan(coef_plan);
		fftw_free(ccoef);
		fftw_free(dcoef);
	}
    double* dcoef; //[norder2]
	fftwf_complex *ccoef;
	fftwf_plan coef_plan;
};
#endif

}
#include "math_chebyshev_def.h"

#endif
