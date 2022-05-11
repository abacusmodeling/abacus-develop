#ifndef STO_CHEBYCHEV_H
#define STO_CHEBYCHEV_H
#include <complex>
#include "fftw3.h"

//----------------------------------------------
// Chebychev Filtering
//----------------------------------------------

class Stochastic_Chebychev
{

	public:

    // constructor and deconstructor
    Stochastic_Chebychev();
    ~Stochastic_Chebychev();

    void init(int ,int);
    
    void calcoef(double fun(double));
	void calcoefc(double fun1(double),double fun2(double));

    double sumallterms();
    
    void calfinalvec(
		void fun(std::complex<double> *in, std::complex<double> *out, const int), 
		std::complex<double> *wavein, 
		std::complex<double> *waveout, 
		const int m = 1);
	void calfinalvec2(
		void fun(std::complex<double> *in, std::complex<double> *out, const int), 
		std::complex<double> *wavein, 
		std::complex<double> *waveout, 
		const int m = 1);

	void calfinalvec_complex(
		void fun(std::complex<double> *in, std::complex<double> *out, const int), 
		std::complex<double> *wavein, 
		std::complex<double> *waveout, 
		const int m = 1);

    bool checkconverge(
		void tfun(std::complex<double> *in, std::complex<double> *out, const int), 
 		std::complex<double> *wavein,
		double& tmax, 
		double& tmin, 
		double stept);

    void calpolyval(
		void fun(std::complex<double> *in, std::complex<double> *out, const int), 
		std::complex<double> *wavein, 
		const int m =1);
	void calpolyvec(
		void fun(std::complex<double> *in, std::complex<double> *out, const int), 
		std::complex<double> *wavein, std::complex<double> *polyvec,
		const int m =1);
	void calpolyval(double x);
	

    int norder;
    int extend;
    int norder2;  // 2 * norder

    double* coef; 
	std::complex<double> * coefc; 

	// expansion coefficient of each order
	// only first norder coefficients are usefull
    double* dcoef; //[norder2]

	// temporary std::complex expansion coefficient of each order
	// only first norder coefficients are usefull
    fftw_complex *ccoef;  //[norder2]

    double *polyvalue;

    fftw_plan plancoef;

    bool initplan;
	bool initcoef;
	bool getcoef;
	bool getpolyval;
	int ndmin; //dim of vector
	int ndmax; //the distance between two closed vectors

    private:

    double ddot_real(const int &m,
    const std::complex<double>* psi_L,
    const std::complex<double>* psi_R);

    template<class T>
    void recurs(
		T *arraynp1, 
		T* arrayn, 
		T *arrayn_1, 
		void fun(T *in,T *out, const int),
		const int m);

};


template<class T>
void Stochastic_Chebychev::recurs(
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
        	arraynp1[i+ib*ndmax]=2*arraynp1[i+ib*ndmax]-arrayn_1[i+ib*ndmax];
    	}
	}
}
#endif// Eelectrons_Chebychev
