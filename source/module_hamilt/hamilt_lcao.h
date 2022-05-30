#ifndef HAMILTLCAO_H
#define HAMILTLCAO_H

#include "hamilt.h"
#include "src_lcao/LCAO_gen_fixedH.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/gint_gamma.h"
#include "src_lcao/gint_k.h"

namespace hamilt
{

// template first for type of k space H matrix elements
// template second for type of temporary matrix, gamma_only fix-gamma-matrix + S-gamma, multi-k fix-Real + S-Real
template <typename T> class HamiltLCAO : public Hamilt
{
  public:
    HamiltLCAO(Gint_Gamma* GG_in, LCAO_gen_fixedH* genH_in, LCAO_Matrix* LM_in)
    {
        this->GG = GG_in;
        this->genH = genH_in;
        this->LM = LM_in;
        this->classname = "HamiltLCAO";
    }
    HamiltLCAO(Gint_k* GK_in, LCAO_gen_fixedH* genH_in, LCAO_Matrix* LM_in)
    {
        this->GK = GK_in;
        this->genH = genH_in;
        this->LM = LM_in;
        this->classname = "HamiltLCAO";
    }
    //~HamiltLCAO();

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const psi::Psi<std::complex<double>>& psi, psi::Psi<std::complex<double>>& hpsi) const override
    {
        // should be updated for iterative diagonalization method
        return;
    };

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    //void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in) override;
    void matrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in) override;

  private:
    void ch_mock();
    void hk_fixed_mock(const int ik);
    void hk_update_mock(const int ik);

    // specific code for matrix()
    void getMatrix(MatrixBlock<T>& hk_in, MatrixBlock<T>& sk_in);

    // there are H and S matrix for each k point in reciprocal space
    // type double for gamma_only case, type complex<double> for multi-k-points case
    std::vector<T> hmatrix_k;
    std::vector<T> smatrix_k;

    void constructFixedReal();
    void constructUpdateReal();

    // temporary class members
    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH* genH = nullptr;

    LCAO_Matrix* LM = nullptr;
};

} // namespace hamilt

#endif