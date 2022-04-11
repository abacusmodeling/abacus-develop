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
void HamiltLCAO<T,T1>::hk_mock()
{
    this->kM.hloc[0] = (T)0.0;
    
    this->kM.sloc[0] = (T)0.0;
}

// complex<double> case for matrix
void HamiltLCAO<std::complex<double>, std::complex<double>>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}
void HamiltLCAO<std::complex<double>, double>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

// double case for matrix
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