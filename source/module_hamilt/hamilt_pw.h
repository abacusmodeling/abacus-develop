#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "hamilt.h"

namespace hamilt
{

class HamiltPW : public Hamilt
{
  public:
    HamiltPW()
    {
      this->classname = "HamiltPW";
    }

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const override;
    virtual void sPsi(const std::complex<double> *psi_in, std::complex<double> *spsi, const size_t size) const override;



  private:

    int current_ik=0;

    int current_spin=0;

    int current_npw=0;

    int max_npw=0;

    const double* current_veff = nullptr;

    void add_nonlocal_pp(
      std::complex<double> *hpsi_in,
      const std::complex<double> *becp,
      const int m) const;
};

} // namespace hamilt

#endif