#include "hamilt_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"

namespace hamilt
{
// case for nspin<4, gamma-only k-point
template class HamiltLCAO<double>;
// case for nspin<4, multi-k-points
// case for nspin == 4, non-collinear spin case
template class HamiltLCAO<std::complex<double>>;

template <typename T> void HamiltLCAO<T>::hk_fixed_mock(const int ik)
{

    // folding_fixedH() should be refactored to there,
    // but now deepks code in this function, should be moved to another place
    return;
}

template <typename T> void HamiltLCAO<T>::hk_update_mock(const int ik)
{
    // update_Hk and update_Hgamma should be refactored to there
    return;
}

template <typename T> void HamiltLCAO<T>::getMatrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in)
{
    hk_in = MatrixBlock<T>{hmatrix_k.data(),
                           (size_t)this->LM->ParaV->nrow,
                           (size_t)this->LM->ParaV->ncol,
                           this->LM->ParaV->desc};
    sk_in = MatrixBlock<T>{smatrix_k.data(),
                           (size_t)this->LM->ParaV->nrow,
                           (size_t)this->LM->ParaV->ncol,
                           this->LM->ParaV->desc};
}

// case for nspin==4
/*template <>
void HamiltLCAO<std::complex<double>, std::complex<double>>::matrix(MatrixBlock<std::complex<double>> &hk_in,
                                                      MatrixBlock<std::complex<double>> &sk_in)
{
    this->getMatrix(hk_in, sk_in);
}*/
// case for nspin<4, multi-k-points
template <>
void HamiltLCAO<std::complex<double>>::matrix(MatrixBlock<std::complex<double>> &hk_in,
                                                      MatrixBlock<std::complex<double>> &sk_in)
{
    this->getMatrix(hk_in, sk_in);
}

// case for nspin<4, gamma_only
template <> void HamiltLCAO<double>::matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in)
{
    this->getMatrix(hk_in, sk_in);
}

template <typename T> void HamiltLCAO<T>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
}

} // namespace hamilt