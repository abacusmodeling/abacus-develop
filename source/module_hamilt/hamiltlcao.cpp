#include "hamiltlcao.h"

namespace ModuleHamilt
{

template<typename T, typename T1>
void HamiltLCAO<T,T1>::ch_mock()
{
    this->kM.resize(1);
    return;
}

template<typename T, typename T1>
void HamiltLCAO<T,T1>::hk_fixed_mock(const int ik)
{
	//folding_fixedH() should be refactored to there, but now deepks code in it
	return;
}

template<typename T, typename T1>
void HamiltLCAO<T, T1>::hk_update_mock(const int ik)
{
    //update_Hk and update_Hgamma should be refactored to there
    return;
}

// case for nspin==4
void HamiltLCAO<std::complex<double>, std::complex<double>>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}
// case for nspin<4, multi-k-points 
void HamiltLCAO<std::complex<double>, double>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

// case for nspin<4, gamma_only
void HamiltLCAO<double, double>::matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

template<typename T, typename T1>
void HamiltLCAO<T, T1>::getMatrix(MatrixBlock<T> hk_in, MatrixBlock<T> sk_in)const
{
    hk_in.p = (T*)kM.hloc.data();
    sk_in.p = (T*)kM.sloc.data();
}

}