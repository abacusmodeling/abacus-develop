#ifndef ELECSTATEPW_H
#define ELECSTATEPW_H

#include "elecstate.h"
#include "module_pw/pw_basis_k.h"

namespace elecstate
{

class ElecStatePW : public ElecState
{
  public:
    ElecStatePW(ModulePW::PW_Basis_K *wfc_basis_in, Charge* chg_in, K_Vectors *pkv_in, int nbands_in) : basis(wfc_basis_in)
    {
        init(chg_in, pkv_in, pkv_in->nks, nbands_in);
        this->classname = "ElecStatePW";
    }
    // void init(Charge* chg_in):charge(chg_in){} override;

    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<std::complex<double>>& psi) override;
    // return current electronic density rho, as a input for constructing Hamiltonian
    // const double* getRho(int spin) const override;

    // update charge density for next scf step
    // void getNewRho() override;

  protected:
    ModulePW::PW_Basis_K *basis;

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    void updateRhoK(const psi::Psi<std::complex<double>>& psi); // override;
    // sum over all pools for rho and ebands
    void parallelK();
    // calcualte rho for each k
    void rhoBandK(const psi::Psi<std::complex<double>>& psi);
};

} // namespace elecstate

#endif