#ifndef HAMILT_PW_HIP_H
#define HAMILT_PW_HIP_H

#include "tools.h"

#ifdef __ROCM
#include "hipfft.h"
#include "hipblas.h"
// #include "rocsolver.h"
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
        ModuleBase::ComplexMatrix &psi,
        const int &nband,
        double *em,
        double *err
    );

	void init_k(const int ik);

	friend class Diago_David;
	// friend class Diago_CG;
    template<class T, class T2, class T3> friend class Diago_CG_CUDA;
	friend class Exx_Lip;
	friend class Hamilt;
    friend class Stochastic_Iter;

	void diagH_subspace(const int ik,
                  const int nstart,
                  const int nbnd,
                  const ModuleBase::ComplexMatrix &psi,
                  ModuleBase::ComplexMatrix &evc,
                  double *en);
    
    void diagH_subspace_cuda(
        const int ik,
        const int nstart,
        const int n_band,
        const hipblasDoubleComplex *psi_c, // matrix
        hipblasDoubleComplex *evc, // matrix
        double *en,
        hipblasDoubleComplex *vkb_c);

    void h_1psi_cuda(
        const int npw,
        const hipblasComplex *psi1d,
        hipblasComplex *hpsi,
        hipblasComplex *spsi,
        hipblasComplex *vkb_c);
    
    void h_1psi_cuda(
        const int npw,
        const hipblasDoubleComplex *psi1d,
        hipblasDoubleComplex *hpsi,
        hipblasDoubleComplex *spsi,
        hipblasDoubleComplex *vkb_c);

    void h_1psi(
        const int npw,
        const std::complex<double> *psi1d,
        std::complex<double> *hpsi,
        std::complex<double> *spsi);

    void h_psi_cuda(
		const hipblasComplex *psi,
		hipblasComplex *hpsi,
        hipblasComplex *vkb_c,
		const int m = 1); 
    
    void h_psi_cuda(
		const hipblasDoubleComplex *psi,
		hipblasDoubleComplex *hpsi,
        hipblasDoubleComplex *vkb_c,
		const int m = 1); 

    void s_1psi_cuda(
        const int npw,
        const hipblasComplex *psi,
        hipblasComplex *spsi);
    
    void s_1psi_cuda(
        const int npw,
        const hipblasDoubleComplex *psi,
        hipblasDoubleComplex *spsi);

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

#ifdef __ROCM
    int *GR_index_d;
    hipblasHandle_t hpw_handle;
#endif

	// add contributions of h*psi from
	// non-local pseduopotentials
	void add_nonlocal_pp_cuda(
		hipblasComplex *hpsi_in,
		const hipblasComplex *becp,
		const hipblasComplex *d_vkb_c,
		const int m);
    
    void add_nonlocal_pp_cuda(
		hipblasDoubleComplex *hpsi_in,
		const hipblasDoubleComplex *becp,
		const hipblasDoubleComplex *d_vkb_c,
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
                          const ModuleBase::ComplexMatrix &psi,
                          const int & m,
                          const std::complex<double> *psik )const ;

	private:

    void diag_zheev
    (
        const int& npw,
        ModuleBase::ComplexMatrix& psi,
        const int& nband,
        double *em,
        double *err
    ) ;

};

#endif
