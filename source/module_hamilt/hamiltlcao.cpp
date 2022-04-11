#include "hamiltlcao.h"

namespace ModuleHamilt
{

template<typename T>
void HamiltLCAO<T>::ch_mock()
{
    this->localM.resize(1);
    return;
}

template<typename T>
void HamiltLCAO<T>::hk_mock()
{
    this->localM.hloc[0] = (T)0.0;
    
    this->localM.sloc[0] = (T)0.0;
}

// complex<double> case for matrix
void HamiltLCAO<std::complex<double>>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

// double case for matrix
void HamiltLCAO<double>::matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

template<typename T>
void HamiltLCAO<T>::getMatrix(MatrixBlock<T> hk_in, MatrixBlock<T> sk_in)const
{
    hk_in.p = (T*)localM.hloc.data();
    sk_in.p = (T*)localM.sloc.data();
}

}