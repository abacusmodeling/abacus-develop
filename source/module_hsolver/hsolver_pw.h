#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "src_pw/pw_basis.h"

namespace hsolver
{

class HSolverPW : public HSolver
{
  public:
    HSolverPW(const PW_Basis* pbas_in)
    {
        this->pbas = pbas_in;
        this->classname = "HSolverPW";
        /*this->init(pbas_in);*/
    }

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/

    void solve(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

  private:
    void hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue);

    void updatePsiK(psi::Psi<std::complex<double>>& psi, const int ik);

    const PW_Basis* pbas = nullptr;

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<double> &h_diag, const int npw, const double* g2kin);

    std::vector<double> precondition;
};

} // namespace hsolver

#endif