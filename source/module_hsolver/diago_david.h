//==========================================================
// AUTHOR : wangjp
// Data :2009-04
// Last Update:
//
// 09-05-10 modify SchmitOrth() diag_zhegvx() as static
// member function
//==========================================================

#ifndef DIAGODAVID_H
#define DIAGODAVID_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "src_pw/structure_factor.h"

namespace hsolver
{

class DiagoDavid : public DiagH
{
  public:
    DiagoDavid(const double* precondition_in);

    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt* phm_in, psi::Psi<std::complex<double>>& phi, double* eigenvalue_in) override;

    static int PW_DIAG_NDIM;

  private:
    int test_david = 0;

    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    void cal_grad(hamilt::Hamilt* phm_in,
                  const int& npw,
                  const int& nbase,
                  const int& notconv,
                  psi::Psi<std::complex<double>>& basis,
                  ModuleBase::ComplexMatrix& hp,
                  ModuleBase::ComplexMatrix& sp,
                  const ModuleBase::ComplexMatrix& vc,
                  const int* unconv,
                  const double* en);

    void cal_elem(const int& npw,
                  int& nbase,
                  const int& notconv,
                  const psi::Psi<std::complex<double>>& basis,
                  const ModuleBase::ComplexMatrix& hp,
                  const ModuleBase::ComplexMatrix& sp,
                  ModuleBase::ComplexMatrix& hc,
                  ModuleBase::ComplexMatrix& sc);

    void refresh(const int& npw,
                 const int& nband,
                 int& nbase,
                 const double* en,
                 const psi::Psi<std::complex<double>>& psi,
                 psi::Psi<std::complex<double>>& basis,
                 ModuleBase::ComplexMatrix& hp,
                 ModuleBase::ComplexMatrix& sp,
                 ModuleBase::ComplexMatrix& hc,
                 ModuleBase::ComplexMatrix& sc,
                 ModuleBase::ComplexMatrix& vc);

    void cal_err(const int& npw,
                 const int& nband,
                 const int& nbase,
                 const ModuleBase::ComplexMatrix& vc,
                 const ModuleBase::ComplexMatrix& hp,
                 const psi::Psi<std::complex<double>>& basis,
                 const double* en,
                 std::complex<double>* respsi);

    void SchmitOrth(const int& npw,
                    const int n_band,
                    const int m,
                    psi::Psi<std::complex<double>>& psi,
                    const ModuleBase::ComplexMatrix& spsi,
                    std::complex<double>* lagrange_m,
                    const int mm_size,
                    const int mv_size);
    void planSchmitOrth(
                    const int nband,
                    int* pre_matrix_mm_m,
                    int* pre_matrix_mv_m);

    void diag_zhegvx(const int& n,
                     const int& m,
                     const ModuleBase::ComplexMatrix& hc,
                     const ModuleBase::ComplexMatrix& sc,
                     const int& ldh,
                     double* e,
                     ModuleBase::ComplexMatrix& vc);

    void diag_mock(hamilt::Hamilt* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in);

    const double* precondition = nullptr;
};

} // namespace hsolver

#endif
