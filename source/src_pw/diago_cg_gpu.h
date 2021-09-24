// #include "ft_gpu.h"
// #include "tools.h"
# ifndef DIAGO_CG_GPU_H
# define DIAGO_CG_GPU_H
#include "cufft.h"
#include "cublas_v2.h"

template<class T, class T2>
class Diago_CG_CUDA
{
public:

    Diago_CG_CUDA();
    ~Diago_CG_CUDA();

    static int moved;

    float ddot_real(
        const int & dim,
        const float2* psi_L,
        const float2* psi_R,
        const bool reduce = true) ;

    double ddot_real(
        const int & dim,
        const double2* psi_L,
        const double2* psi_R,
        const bool reduce = true) ;

    float2 ddot(
        const int & dim,
        const float2* psi_L,
        const float2* psi_R ) ;

    double2 ddot(
        const int & dim,
        const double2* psi_L,
        const double2* psi_R ) ;

    float2 ddot(
        const int & dim,
        const float2 *psi, // matrix
        const int & m,
        float2 *psik ) ;

    double2 ddot(
        const int & dim,
        const double2 *psi, // matrix
        const int & m,
        double2 *psik ) ;

    float2 ddot(
        const int & dim,
        const float2 *psi_L, // matrix
        const int & m,
        const float2 *psi_R, // matrix
        const int & n) ;
    
    double2 ddot(
        const int & dim,
        const double2 *psi_L, // matrix
        const int & m,
        const double2 *psi_R, // matrix
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
        double &avg_iter);

    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const float2 *psi, // matrix
        float2 *spsi,
        float2 *psi_m
    );

    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const double2 *psi, // matrix
        double2 *spsi,
        double2 *psi_m
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
        float2 *g,
        float2 *sg,
        float2 *lagrange,
        const float2 *eigenfunction, // matrix
        const int m);
    
    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        double2 *g,
        double2 *sg,
        double2 *lagrange,
        const double2 *eigenfunction, // matrix
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
