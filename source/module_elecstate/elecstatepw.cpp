#include "elecstatepw.h"

namespace ModuleElecS
{

void ElecStatePW::init(
    Charge_Broyden* chr_in)
{
    this->chr = chr_in;
}

const MatrixBlock<double> ElecStatePW::getRho()const
{
    MatrixBlock<double> temp{&(this->chr->rho[0][0]), 1,1};//this->chr->get_nspin(), this->chr->get_nrxx()};
    return temp;
}

void ElecStatePW::updateRhoK(const ModulePsi::Psi<std::complex<double>> &psi)
{
    //this->chr->sum_band_k(psi);
    return;
}

void ElecStatePW::getNewRho()
{
    return;
}

}