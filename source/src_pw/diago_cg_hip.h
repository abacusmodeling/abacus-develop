# ifndef DIAGO_CG_CUH
# define DIAGO_CG_CUH

#ifdef __ROCM
#include "hipfft.h"
#include "hipblas.h"
// #include "hip/hip_complex.h"
#endif

template<class T, class T2, class T3>
class Diago_CG_CUDA
{
public:

    Diago_CG_CUDA();
    ~Diago_CG_CUDA();

    static int moved;

    float ddot_real(
        const int & dim,
        const hipblasComplex* psi_L,
        const hipblasComplex* psi_R,
        const bool reduce = true) ;

    double ddot_real(
        const int & dim,
        const hipblasDoubleComplex* psi_L,
        const hipblasDoubleComplex* psi_R,
        const bool reduce = true) ;

    hipblasComplex ddot(
        const int & dim,
        const hipblasComplex* psi_L,
        const hipblasComplex* psi_R ) ;

    hipblasDoubleComplex ddot(
        const int & dim,
        const hipblasDoubleComplex* psi_L,
        const hipblasDoubleComplex* psi_R ) ;

    hipblasComplex ddot(
        const int & dim,
        const hipblasComplex *psi, // matrix
        const int & m,
        hipblasComplex *psik ) ;

    hipblasDoubleComplex ddot(
        const int & dim,
        const hipblasDoubleComplex *psi, // matrix
        const int & m,
        hipblasDoubleComplex *psik ) ;

    hipblasComplex ddot(
        const int & dim,
        const hipblasComplex *psi_L, // matrix
        const int & m,
        const hipblasComplex *psi_R, // matrix
        const int & n) ;
    
    hipblasDoubleComplex ddot(
        const int & dim,
        const hipblasDoubleComplex *psi_L, // matrix
        const int & m,
        const hipblasDoubleComplex *psi_R, // matrix
        const int & n) ;

    void diag(
        T2 *phi, // matrix
        T *e,
        T2 *vkb_c,
        const int &dim,
        const int &dmx,
        const int &n_band,
        const T *precondition,
        const T &eps,
        const int &maxter,
        const bool &reorder,
        int &notconv,
        double &avg_iter);

    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const hipblasComplex *psi, // matrix
        hipblasComplex *spsi,
        hipblasComplex *psi_m
    );

    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const hipblasDoubleComplex *psi, // matrix
        hipblasDoubleComplex *spsi,
        hipblasDoubleComplex *psi_m
    );

private:

    int test_cg;

#ifdef __ROCM
    hipblasHandle_t diag_handle;
#endif

    void calculate_gradient(
        const T* precondition,
        const int dim,
        const T2 *hpsi,
        const T2 *spsi,
        T2 *g,
        T2 *pspsi);

    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        hipblasComplex *g,
        hipblasComplex *sg,
        hipblasComplex *lagrange,
        const hipblasComplex *eigenfunction, // matrix
        const int m);
    
    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        hipblasDoubleComplex *g,
        hipblasDoubleComplex *sg,
        hipblasDoubleComplex *lagrange,
        const hipblasDoubleComplex *eigenfunction, // matrix
        const int m);

    void calculate_gamma_cg(
        const int iter,
        const int dim,
        const T *precondition,
        const T2 *g,
        const T2 *sg,
        T2 *psg,
        T2 *cg,
        T &gg_last,
        const T &cg0,
        const T &theta,
        const T2 *psi_m);

    bool update_psi(
        const int dim,
        T &cg_norm,
        T &theta,
        T2 *hcg,
        const T2 *cg,
        T2 *scg,
        T2 *psi_m ,
        T &eigenvalue,
        const T &threshold,
        T2 *hpsi,
        T2 *spsi,
        T2 *vkb_c);

};
# endif
