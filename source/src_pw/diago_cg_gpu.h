// #include "ft_gpu.h"
// #include "tools.h"
# ifndef DIAGO_CG_GPU_H
# define DIAGO_CG_GPU_H
#include "cufft.h"
#include "cublas_v2.h"

typedef cufftDoubleComplex CUFFT_COMPLEX;

class Diago_CG_GPU
{
public:

    Diago_CG_GPU();
    ~Diago_CG_GPU();

    static int moved;

    double ddot_real(
        const int & dim,
        const CUFFT_COMPLEX* psi_L,
        const CUFFT_COMPLEX* psi_R,
        const bool reduce = true) ;

    CUFFT_COMPLEX ddot(
        const int & dim,
        const CUFFT_COMPLEX* psi_L,
        const CUFFT_COMPLEX* psi_R ) ;

    CUFFT_COMPLEX ddot(
        const int & dim,
        const CUFFT_COMPLEX *psi, // matrix
        const int & m,
        CUFFT_COMPLEX *psik ) ;

    CUFFT_COMPLEX ddot(
        const int & dim,
        const CUFFT_COMPLEX *psi_L, // matrix
        const int & m,
        const CUFFT_COMPLEX *psi_R, // matrix
        const int & n) ;

    void diag(
        CUFFT_COMPLEX *phi, // matrix
        double *e,
        const int &dim,
        const int &dmx,
        const int &n_band,
        const double *precondition,
        const double &eps,
        const int &maxter,
        const bool &reorder,
        int &notconv,
        double &avg_iter);

    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const CUFFT_COMPLEX *psi, // matrix
        CUFFT_COMPLEX *spsi,
        CUFFT_COMPLEX *psi_m
    );

private:

    int test_cg;
    cublasHandle_t diag_handle;

    void calculate_gradient(
        const double* precondition,
        const int dim,
        const CUFFT_COMPLEX *hpsi,
        const CUFFT_COMPLEX *spsi,
        CUFFT_COMPLEX *g,
        CUFFT_COMPLEX *pspsi);

    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        CUFFT_COMPLEX *g,
        CUFFT_COMPLEX *sg,
        CUFFT_COMPLEX *lagrange,
        const CUFFT_COMPLEX *eigenfunction, // matrix
        const int m);

    void calculate_gamma_cg(
        const int iter,
        const int dim,
        const double *precondition,
        const CUFFT_COMPLEX *g,
        const CUFFT_COMPLEX *sg,
        CUFFT_COMPLEX *psg,
        CUFFT_COMPLEX *cg,
        double &gg_last,
        const double &cg0,
        const double &theta,
        const CUFFT_COMPLEX *psi_m);

    bool update_psi(
        const int dim,
        double &cg_norm,
        double &theta,
        CUFFT_COMPLEX *hcg,
        const CUFFT_COMPLEX *cg,
        CUFFT_COMPLEX *scg,
        CUFFT_COMPLEX *psi_m ,
        double &eigenvalue,
        const double &threshold,
        CUFFT_COMPLEX *hpsi,
        CUFFT_COMPLEX *spsi);

};
# endif
