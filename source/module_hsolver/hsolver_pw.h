#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "module_pw/pw_basis_k.h"

namespace hsolver
{

class HSolverPW : public HSolver
{
  public:
    HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in)
    {
        this->wfc_basis = wfc_basis_in;
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

    ModulePW::PW_Basis_K* wfc_basis = nullptr;

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<double> &h_diag, const int ik, const int npw);

    std::vector<double> precondition;
};

} // namespace hsolver

#endif