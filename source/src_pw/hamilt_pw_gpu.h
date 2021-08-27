#ifndef HAMILT_PW_GPU_H
#define HAMILT_PW_GPU_H

#include "tools.h"

#ifdef __CUDA
#include "cufft.h"
#include "cublas_v2.h"
typedef cufftDoubleComplex CUFFT_COMPLEX;
#endif

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

	friend class Diago_David;
	// friend class Diago_CG;
    friend class Diago_CG_GPU;
	friend class Exx_Lip;
	friend class Hamilt;
    friend class Stochastic_Iter;

	void diagH_subspace(const int ik,
                  const int nstart,
                  const int nbnd,
                  const ComplexMatrix &psi,
                  ComplexMatrix &evc,
                  double *en);

    void h_1psi_gpu(
        const int npw,
        const CUFFT_COMPLEX *psi1d,
        CUFFT_COMPLEX *hpsi,
        CUFFT_COMPLEX *spsi);

    void h_1psi(
        const int npw,
        const std::complex<double> *psi1d,
        std::complex<double> *hpsi,
        std::complex<double> *spsi);

    void h_psi_gpu(
		const CUFFT_COMPLEX *psi,
		CUFFT_COMPLEX *hpsi,
		const int m = 1); // qianrui add a default parameter 2021-3-31

    void s_1psi_gpu(
        const int npw,
        const CUFFT_COMPLEX *psi,
        CUFFT_COMPLEX *spsi);

    void s_1psi(
        const int npw,
        const std::complex < double> *psi,
        std::complex < double> *spsi);

    void h_psi(
		const std::complex<double> *psi,
		std::complex<double> *hpsi,
		const int m = 1);

	private:

    int *GR_index;

#ifdef __CUDA
    int *GR_index_d;
#endif

	// add contributions of h*psi from
	// non-local pseduopotentials
	void add_nonlocal_pp_gpu(
		CUFFT_COMPLEX *hpsi_in,
		const CUFFT_COMPLEX *becp,
		const CUFFT_COMPLEX *d_vkb_c,
		const int m);

    void add_nonlocal_pp(
		std::complex<double> *hpsi,
		const std::complex<double> *becp,
		const int m);

	private:

	double ddot_real(
		const int& npw,
		const std::complex<double>* psi_L,
		const std::complex<double>* psi_R)const;

    std::complex<double> ddot( const int& npw,
                          const std::complex<double> * psi_L,
                          const std::complex<double> * psi_R )const ;

    std::complex<double> just_ddot( const int& npw,
                          const std::complex<double> * psi_L,
                          const std::complex<double> * psi_R )const ;

    std::complex<double> ddot( const int & npw,
                          const ComplexMatrix &psi,
                          const int & m,
                          const std::complex<double> *psik )const ;

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
