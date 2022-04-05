#ifndef DIAGCG_H
#define DIAGCG_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_psi/psi.h"

#include "src_pw/hamilt_pw.h"
#include "src_pw/pw_basis.h"

namespace ModuleHSolver
{

class DiagoCG : public DiagH
{
	public:

    DiagoCG(
        Hamilt_PW* hpw_in, 
        const PW_Basis* pbas_in, 
        const double *precondition_in);
    ~DiagoCG();

    //virtual void init(){};
    virtual int diag (
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<std::complex<double>> &phi,
        double *eigenvalue_in) override
    {
        return this->diag_mock(phi, eigenvalue_in);
    }

	private:
    /// static variables, used for passing control variables 
    /// max steps for cg diagonalization iteration, read in from INPUT file.
    static int maxter;
    /// if eigenvalue and eigenvectors should be reordered after diagonalization, it is always be true.
    static bool reorder;
    /// threshold used in cg, read int from INPUT file.
    static double eps;
    /// average steps of last cg diagonalization for each band.
    static double avg_iter;

    /// temp operator pointer
    Hamilt_PW* hpw;
    const PW_Basis* pbas;

    int test_cg;

    /// inside variables and vectors, used by inside functions.
    /// row size for input psi matrix
    int n_band;
    /// col size for input psi matrix
    int dmx;
    /// non-zero col size for inputted psi matrix 
    int dim;
    /// precondition for cg diag
    const double* precondition;
    /// eigenvalue results
    double* eigenvalue;
    

    /// temp vector for S|psi> for one band, size dim
    std::vector<std::complex<double>> sphi;
    /// temp vector for , size dim
    std::vector<std::complex<double>> scg;
    /// temp vector for H|psi> for one band, size dim
    std::vector<std::complex<double>> hphi; 
    /// temp vector for , size dim
    std::vector<std::complex<double>> gradient;
    /// temp vector for , size dim
    std::vector<std::complex<double>> cg;
    /// temp vector for , size dim   
    std::vector<std::complex<double>> g0;
    /// temp vector for store psi in sorting with eigenvalues, size dim  
    std::vector<std::complex<double>> pphi; 
    /// temp vector for matrix eigenvector * vector S|psi> , size m_band
    std::vector<std::complex<double>> lagrange; 
    /// temp vector for new psi for one band, size dim
    std::vector<std::complex<double>> phi_m;

    void calculate_gradient();

    void orthogonal_gradient(
        const ModulePsi::Psi<std::complex<double>> &eigenfunction,
        const int m);

    void calculate_gamma_cg(
        const int iter,
        double &gg_last,
        const double &cg0,
        const double &theta);

    bool update_psi(
        double &cg_norm,
        double &theta,
        double &eigenvalue);

    void schmit_orth(
        const int &m,
        const ModulePsi::Psi<std::complex<double>> &psi
    );

    //used in diag() for template replace Hamilt with Hamilt_PW
    int diag_mock(ModulePsi::Psi<std::complex<double>> &phi,
        double *eigenvalue_in);

};

}
#endif
