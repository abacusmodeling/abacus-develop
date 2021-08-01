#ifndef HAMILT_PW_GPU_H
#define HAMILT_PW_GPU_H

#include "tools.h"
#include "cufft.h"
#include "cublas_v2.h"

typedef cufftDoubleComplex CUFFT_COMPLEX;

class Hamilt_PW_GPU
{

public:

    Hamilt_PW_GPU();
    ~Hamilt_PW_GPU();

    static int moved;

	void allocate(
		const int &npwx, // number of plane wave (max)
		const int &npol, // polarization 
		const int &nkb,  // number of non-local pseudopotential projectors 
		const int &nrxx); // number of grids on this processor

    void init_k(const int ik);

	private:
	
	friend class Diago_David;
	friend class Diago_CG;
    friend class Diago_CG_GPU;
	friend class Exx_Lip;
	friend class Hamilt;
    friend class Stochastic_Iter;

    void h_1psi(
        const int npw,
        const CUFFT_COMPLEX *psi1d,
        CUFFT_COMPLEX *hpsi,
        CUFFT_COMPLEX *spsi);

    void h_psi( 
		const CUFFT_COMPLEX *psi, 
		CUFFT_COMPLEX *hpsi, 
		const int m = 1); // qianrui add a default parameter 2021-3-31

    void s_1psi(
        const int npw,
        const CUFFT_COMPLEX *psi,
        CUFFT_COMPLEX *spsi);

	private:

    int *GR_index;

	// add contributions of h*psi from
	// non-local pseduopotentials
    void add_nonlocal_pp(
		complex<double> *hpsi, 
		const complex<double> *becp, 
		const int m);

	private:

};

#endif 
