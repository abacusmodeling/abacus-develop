// #include "ft_gpu.h"
// #include "tools.h"
# ifndef DIAGO_CG_GPU_H
# define DIAGO_CG_GPU_H
#include "cufft.h"
#include "cublas_v2.h"

// typedef cufftDoubleComplex CUFFT_COMPLEX;

template<class T, class T2>
class Diago_CG_GPU
{
public:

    Diago_CG_GPU();
    ~Diago_CG_GPU();

    static int moved;

    T ddot_real(
        const int & dim,
        const T2* psi_L,
        const T2* psi_R,
        const bool reduce = true) ;

    T2 ddot(
        const int & dim,
        const T2* psi_L,
        const T2* psi_R ) ;

    T2 ddot(
        const int & dim,
        const T2 *psi, // matrix
        const int & m,
        T2 *psik ) ;

    T2 ddot(
        const int & dim,
        const T2 *psi_L, // matrix
        const int & m,
        const T2 *psi_R, // matrix
        const int & n) ;

    void diag(
        T2 *phi, // matrix
        T *e,
        const int &dim,
        const int &dmx,
        const int &n_band,
        const T *precondition,
        const T &eps,
        const int &maxter,
        const bool &reorder,
        int &notconv,
        T &avg_iter);

    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const T2 *psi, // matrix
        T2 *spsi,
        T2 *psi_m
    );

private:

    int test_cg;
    cublasHandle_t diag_handle;

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
        T2 *g,
        T2 *sg,
        T2 *lagrange,
        const T2 *eigenfunction, // matrix
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
        T2 *spsi);

};
# endif
