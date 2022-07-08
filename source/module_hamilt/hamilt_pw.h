#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "hamilt.h"
#if ((defined __CUDA) || (defined __ROCM))

#ifdef __CUDA
#include "src_pw/hamilt_pw.cuh"
#else
#include "src_pw/hamilt_pw_hip.h"
#endif

#else
#include "src_pw/hamilt_pw.h"
#endif

namespace hamilt
{

class HamiltPW : public Hamilt
{
  public:
    HamiltPW(Hamilt_PW* hpw_in)
    {
      this->hpw = hpw_in;
      this->classname = "HamiltPW";
    }

    // construct Hamiltonian matrix with inputed electonic density
    void constructHamilt() override
    {
        this->ch_mock();
    }

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override
    {
        this->hk_mock(ik);
    };

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const psi::Psi<std::complex<double>>& psi, psi::Psi<std::complex<double>>& hpsi) const override
    {
        this->hpsi_mock(psi, hpsi);
    };



  private:
    Hamilt_PW* hpw;
    void ch_mock();
    void hk_mock(const int ik);
    void hpsi_mock(const psi::Psi<std::complex<double>>& psi, psi::Psi<std::complex<double>>& hpsi) const;
};

} // namespace hamilt

#endif