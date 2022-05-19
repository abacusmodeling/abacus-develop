#include "./constants.h"

namespace ModuleBase
{

template<typename REAL>
template<class T>
void Chebyshev<REAL>:: calcoef_real(T *ptr, REAL (T::*fun)(REAL))
{
    complex<REAL> *pcoef = (complex<REAL> *)this->fftw.ccoef;

    //three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)((REAL) cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }

    //this->fftw.dcoef --FFT--> pcoef
    fftw_execute(this->fftw.coef_plan);

    for(int i = 0; i<norder; ++i)
    {
        REAL phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_real[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3;
        }
        else
        {
            coef_real[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3;
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos(i*ModuleBase::TWO_PI/norder2));
    }

    //this->fftw.dcoef --FFT--> pcoef
    fftw_execute(this->fftw.coef_plan);

    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_real[i] += real(pcoef[i]) / norder2 * 1 / 3;
        }
        else
        {
            coef_real[i] += real(pcoef[i]) / norder2 * 2 / 3;
        }
    } 

    getcoef = true;
	return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calcoef_complex(T *ptr, std::complex<REAL> (T::*fun)(std::complex<REAL>))
{
    complex<double> *pcoef = (complex<double> *)this->fftw.ccoef;

    //three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos((i+0.5)*ModuleBase::TWO_PI/norder2)).real();
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos((i+0.5)*ModuleBase::TWO_PI/norder2)).imag();
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos(i*ModuleBase::TWO_PI/norder2)).real();
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos(i*ModuleBase::TWO_PI/norder2)).imag();
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    getcoef = true;
	return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calcoef_pair(T *ptr, REAL (T::*fun1)(REAL), REAL (T::*fun2)(REAL))
{
    complex<double> *pcoef = (complex<double> *)this->fftw.ccoef;

    //three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun1)(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun2)(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        double phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun1)(cos(i*ModuleBase::TWO_PI/norder2));
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun2)(cos(i*ModuleBase::TWO_PI/norder2));
    }
    fftw_execute(this->fftw.coef_plan);
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    getcoef = true;
	return;
}



}