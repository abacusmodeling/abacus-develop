#ifndef DIAGO_CG_H
#define DIAGO_CG_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"

class Diago_CG
{
	public:

    Diago_CG();
    ~Diago_CG();

    static int moved;

    static double ddot_real(
        const int & dim,
        const std::complex<double>* psi_L,
        const std::complex<double>* psi_R,
        const bool reduce = true) ;

    static std::complex<double> ddot(
        const int & dim,
        const std::complex<double>* psi_L,
        const std::complex<double>* psi_R ) ;


    static std::complex<double> ddot(
        const int & dim,
        const ModuleBase::ComplexMatrix &psi,
        const int & m,
        std::complex<double> *psik ) ;

    static std::complex<double> ddot(
        const int & dim,
        const ModuleBase::ComplexMatrix &psi_L,
        const int & m,
        const ModuleBase::ComplexMatrix &psi_R,
        const int & n) ;

    void diag(
        ModuleBase::ComplexMatrix &phi,
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

    static void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &end,
        const ModuleBase::ComplexMatrix &psi,
        std::complex<double> *spsi,
        std::complex<double> *psi_m
    );

	private:

    int test_cg;

    void calculate_gradient(
        const double* precondition,
        const int dim,
        const std::complex<double> *hpsi,
        const std::complex<double> *spsi,
        std::complex<double> *g,
        std::complex<double> *pspsi);

    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        std::complex<double> *g,
        std::complex<double> *scg,
        std::complex<double> *lagrange,
        const ModuleBase::ComplexMatrix &eigenfunction,
        const int m);

    void calculate_gamma_cg(
        const int iter,
        const int dim,
        const double *precondition,
        const std::complex<double> *g,
        const std::complex<double> *sg,
        std::complex<double> *psg,
        std::complex<double> *cg,
        double &gg_last,
        const double &cg0,
        const double &theta,
        const std::complex<double> *psi_m);

    bool update_psi(
        const int dim,
        double &cg_norm,
        double &theta,
        std::complex<double> *hcg,
        const std::complex<double> *cg,
        std::complex<double> *scg,
        std::complex<double> *psi_m ,
        double &eigenvalue,
        const double &threshold,
        std::complex<double> *hpsi,
        std::complex<double> *spsi);

};
#endif
