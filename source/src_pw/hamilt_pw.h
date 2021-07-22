#ifndef HAMILT_PW_H
#define HAMILT_PW_H

#include "tools.h"

class Hamilt_PW
{

public:

    Hamilt_PW();
    ~Hamilt_PW();

    static int moved;

	void allocate(
		const int &npwx, // number of plane wave (max)
		const int &npol, // polarization 
		const int &nkb,  // number of non-local pseudopotential projectors 
		const int &nrxx); // number of grids on this processor

    void cal_err
    (
        const int &npw,
        ComplexMatrix &psi,
        const int &nband,
        double *em,
        double *err
    );

    void init_k(const int ik);

	private:
	
	friend class Diago_David;
	friend class Diago_CG;
	friend class Exx_Lip;
	friend class Hamilt;
    friend class Stochastic_Iter;

    void diagH_subspace(const int ik,
                  const int nstart,
                  const int nbnd,
                  const ComplexMatrix &psi,
                  ComplexMatrix &evc,
                  double *en);

    void h_1psi(
        const int npw,
        const complex<double> *psi1d,
        complex<double> *hpsi,
        complex<double> *spsi);

    void h_psi( 
		const complex<double> *psi, 
		complex<double> *hpsi, 
		const int m = 1); // qianrui add a default parameter 2021-3-31

    void s_1psi(
        const int npw,
        const complex < double> *psi,
        complex < double> *spsi);

	private:

    int *GR_index;

    complex<double> *psi_one;

    // hpsi and spsi
    complex<double> *hpsi;
    complex<double> *spsi;
    complex<double> *Bec;

	// add contributions of h*psi from
	// non-local pseduopotentials
    void add_nonlocal_pp(
		complex<double> *hpsi, 
		const complex<double> *becp, 
		const int m);

	private:

    double ddot_real( 
		const int& npw, 
		const complex<double>* psi_L, 
		const complex<double>* psi_R)const;

    complex<double> ddot( const int& npw,
                          const complex<double> * psi_L,
                          const complex<double> * psi_R )const ;

    complex<double> just_ddot( const int& npw,
                          const complex<double> * psi_L,
                          const complex<double> * psi_R )const ;

    complex<double> ddot( const int & npw,
                          const ComplexMatrix &psi,
                          const int & m,
                          const complex<double> *psik )const ;

	private:

    void diag_zheev
    (
        const int& npw,
        ComplexMatrix& psi,
        const int& nband,
        double *em,
        double *err
    ) ;

};

#endif 
