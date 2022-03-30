#include "psi.h"

namespace ModulePsi
{

template<typename T>
Psi<T>::Psi(
    const int& nk_in, 
    const int& nbd_in, 
    const int& nbs_in, 
    const bool spin_method_in
)
{
    this->resize(nk_in, nbd_in, nbs_in);
    this->spin_method = spin_method_in;
    this->current_b = 0;
    this->current_k = 0;
}

template<typename T>
Psi<T>::Psi(
    const Psi<T>& psi_in, 
    const int& nk_in, 
    const int& nbd_in
)
{
    assert(nk_in<=psi_in.get_nk() && nbd_in<=psi_in.get_nbd());
    this->resize(nk_in, nbd_in, psi_in.get_nbasis());
    //if size of k is 1, copy from Psi in current_k, 
    //else copy from start of Psi
    if(nk_in==1) for(size_t index=0; index<this->size();++index)
    {
        psi[index] = psi_in.get_pointer()[index];
        //current_k for this Psi only keep the spin index same as the copied Psi
        this->current_k = psi_in.get_spin();
    } 
    else for(size_t index=0; index<this->size();++index) psi[index] = psi_in[index];

    this->spin_method = psi_in.spin_method;
}

template<typename T>
void Psi<T>::resize(
    const int nks_in,
    const int nbands_in,
    const int nbasis_in)
{
    assert(nks_in>0 && nbands_in>0 && nbasis_in>0);
    this->psi.resize(nks_in * nbands_in * nbasis_in);
    this->nk = nks_in;
    this->nbands = nbands_in;
    this->nbasis = nbasis_in;
    this->psi_current = psi.data();

    return;
}

}