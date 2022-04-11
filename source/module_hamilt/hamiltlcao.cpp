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
void HamiltLCAO<std::complex<double>>::matrix(const MatrixBlock<std::complex<double>> hk_in, const MatrixBlock<std::complex<double>> sk_in)
{
    hk_in.p[0] = localM.hloc[0];
    sk_in.p[0] = localM.sloc[0];
}

}