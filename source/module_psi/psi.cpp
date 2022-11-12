#include "psi.h"
#include <cassert>
#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
#include "module_psi/include/device.h"

#include <complex>

using namespace psi;

Range::Range(const size_t range_in) 
{
    k_first = 1;
    index_1 = 0;
    range_1 = range_in;
    range_2 = range_in;
}

Range::Range(
    const bool k_first_in, 
    const size_t index_1_in, 
    const size_t range_1_in, 
    const size_t range_2_in)
{
    k_first = k_first_in;
    index_1 = index_1_in;
    range_1 = range_1_in;
    range_2 = range_2_in;
}


template<typename T, typename Device>
Psi<T, Device>::Psi() 
{
    this->npol = GlobalV::NPOL;
    this->device = device::get_device_type<Device>(this->ctx);
}

template<typename T, typename Device>
Psi<T, Device>::~Psi() 
{
    delete_memory_op()(this->ctx, this->psi);
}

template<typename T, typename Device>
Psi<T, Device>::Psi(
    const int* ngk_in) 
{
    this->ngk = ngk_in;
    this->npol = GlobalV::NPOL;
    this->device = device::get_device_type<Device>(this->ctx);
}

template<typename T, typename Device>
Psi<T, Device>::Psi(
    int nk_in, 
    int nbd_in, 
    int nbs_in, 
    const int* ngk_in) 
{
    this->ngk = ngk_in;
    this->current_b = 0;
    this->current_k = 0;
    this->npol = GlobalV::NPOL;
    this->device = device::get_device_type<Device>(this->ctx);
    this->resize(nk_in, nbd_in, nbs_in);
}

template<typename T, typename Device>
Psi<T, Device>::Psi(
    const Psi& psi_in, 
    const int nk_in, 
    int nband_in) 
{
    assert(nk_in<=psi_in.get_nk());
    if(nband_in == 0) {
        nband_in = psi_in.get_nbands();
    }
    this->device = psi_in.device;
    this->resize(nk_in, nband_in, psi_in.get_nbasis());
    this->ngk = psi_in.ngk;
    this->npol = psi_in.npol;    
    if(nband_in <= psi_in.get_nbands()) {
        // copy from Psi from psi_in(current_k, 0, 0), 
        // if size of k is 1, current_k in new Psi is psi_in.current_k 
        if(nk_in == 1) {
            //current_k for this Psi only keep the spin index same as the copied Psi
            this->current_k = psi_in.get_current_k();
        } 
        synchronize_memory_op()(
            this->ctx, 
            psi_in.get_device(),
            this->psi, 
            psi_in.get_pointer(), 
            this->size());
    }
}

template<typename T, typename Device>
Psi<T, Device>::Psi(
    T* psi_pointer, 
    const Psi& psi_in, 
    const int nk_in, 
    int nband_in)
{
    this->device = device::get_device_type<Device>(this->ctx);
    assert(this->device == psi_in.device);
    assert(nk_in<=psi_in.get_nk());
    if(nband_in == 0) {
        nband_in = psi_in.get_nbands();
    }
    this->ngk = psi_in.ngk;
    this->npol = psi_in.npol;
    this->nk = nk_in;
    this->nbands = nband_in;
    this->nbasis = psi_in.nbasis;
    this->psi_current = psi_pointer;
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
Psi<T, Device>::Psi(
    const Psi<T_in, Device_in>& psi_in)
{
    this->ngk = psi_in.get_ngk_pointer();
    this->npol = psi_in.npol;
    this->nk = psi_in.get_nk();
    this->nbands = psi_in.get_nbands();
    this->nbasis = psi_in.get_nbasis();
    this->current_k = psi_in.get_current_k();
    this->current_b = psi_in.get_current_b();
    this->current_nbasis = psi_in.get_current_nbas();
    this->k_first = psi_in.get_k_first();
    // this function will copy psi_in.psi to this->psi no matter the device types of each other.
    this->device = device::get_device_type<Device>(this->ctx);
    this->resize(psi_in.get_nk(), psi_in.get_nbands(), psi_in.get_nbasis());
    memory::synchronize_memory_op<T, Device, Device_in>()(
        this->ctx,
        psi_in.get_device(),
        this->psi,
        psi_in.get_pointer(), 
        psi_in.size());
    this->psi_current = this->psi + psi_in.get_psi_bias();
}

template<typename T, typename Device>
void Psi<T, Device>::resize(
    const int nks_in, 
    const int nbands_in,
    const int nbasis_in)
{
    assert(nks_in>0 && nbands_in>=0 && nbasis_in>0);
    // This function will delete the psi array first(if psi exist), then malloc a new memory for it.
    resize_memory_op()(
        this->ctx,
        this->psi, 
        nks_in * nbands_in * nbasis_in);
    this->nk = nks_in;
    this->nbands = nbands_in;
    this->nbasis = nbasis_in;
    this->current_nbasis = nbasis_in;
    this->psi_current = this->psi;
}

template<typename T, typename Device>
T* Psi<T, Device>::get_pointer()
{
    return this->psi_current;
}

template<typename T, typename Device>
T* Psi<T, Device>::get_pointer(const int& ibands)
{
    assert(ibands >= 0 && ibands < this->nbands);
    return &this->psi_current[ibands * this->nbasis];
}

template<typename T, typename Device>
const T* Psi<T, Device>::get_pointer() const
{
    return this->psi_current;
}

template<typename T, typename Device>
const int* Psi<T, Device>::get_ngk_pointer() const
{
    return this->ngk;
}

template<typename T, typename Device>
const bool& Psi<T, Device>::get_k_first() const
{
    return this->k_first;
}

template<typename T, typename Device>
const Device* Psi<T, Device>::get_device() const
{
    return this->ctx;
}

template<typename T, typename Device>
const int& Psi<T, Device>::get_psi_bias() const
{
    return this->psi_bias;
}

template<typename T, typename Device>
const T* Psi<T, Device>::get_pointer(const int& ibands) const
{
    assert(ibands >= 0 && ibands < this->nbands);
    return &this->psi_current[ibands * this->nbasis];
}

template<typename T, typename Device>
const int& Psi<T, Device>::get_nk() const
{
    return this->nk;
}

template<typename T, typename Device>
const int& Psi<T, Device>::get_nbands() const
{
    return this->nbands;
}

template<typename T, typename Device>
const int& Psi<T, Device>::get_nbasis() const
{
    return this->nbasis;
}

template<typename T, typename Device>
size_t Psi<T, Device>::size() const
{   
    if (this-> psi == nullptr) {
        return 0;
    }
    return this->nk * this->nbands * this->nbasis;
}

template<typename T, typename Device>
void Psi<T, Device>::fix_k(const int ik) const
{
    assert(ik>=0);
    this->current_k = ik;
    if(this->ngk != nullptr && this->npol != 2) this->current_nbasis = this->ngk[ik];
    else this->current_nbasis = this->nbasis;
    this->current_b = 0;
    if(ik >= this->nk) {
        // mem_saver case
        this->psi_current = const_cast<T*>(&(this->psi[0]));
        this->psi_bias = 0;
    }
    else {
        this->psi_current = const_cast<T*>(&(this->psi[ik * this->nbands * this->nbasis]));
        this->psi_bias = ik * this->nbands * this->nbasis;
    }
}

template<typename T, typename Device>
void Psi<T, Device>::fix_band(const int iband) const
{
    assert(iband >= 0 && iband < this->nbands);
    this->current_b = iband;
}

template<typename T, typename Device>
T& Psi<T, Device>::operator()(const int ik, const int ibands, const int ibasis)
{
    assert(ik >= 0 && ik < this->nk);
	assert(ibands >= 0 && ibands < this->nbands);	
    assert(ibasis >= 0 && ibasis < this->nbasis);
    if (this->device != CpuDevice) {
        #if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::WARNING_QUIT("Psi operator ()", "GPU psi cannot fetch value by an overloaded operator ()!");
        #endif 
    }
	return this->psi[(ik * this->nbands + ibands) * this->nbasis + ibasis];
}

template<typename T, typename Device>
const T& Psi<T, Device>::operator()(const int ik, const int ibands, const int ibasis) const
{
    assert(ik >= 0 && ik < this->nk);
	assert(ibands >= 0 && ibands < this->nbands);	
    assert(ibasis >= 0 && ibasis < this->nbasis);
    if (this->device != CpuDevice) {
        #if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::WARNING_QUIT("Psi operator ()", "GPU psi cannot fetch value by an overloaded operator ()!");
        #endif 
    }
	return this->psi[(ik * this->nbands + ibands) * this->nbasis + ibasis];
}

template<typename T, typename Device>
T& Psi<T, Device>::operator()(const int ibands, const int ibasis)
{
    assert(this->current_b == 0);
	assert(ibands >= 0 && ibands < this->nbands);	
    assert(ibasis >= 0 && ibasis < this->nbasis);
    if (this->device != CpuDevice) {
        #if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::WARNING_QUIT("Psi operator ()", "GPU psi cannot fetch value by an overloaded operator ()!");
        #endif 
    }
	return this->psi_current[ibands * this->nbasis + ibasis];
}

template<typename T, typename Device>
const T& Psi<T, Device>::operator()(const int ibands, const int ibasis) const
{
    assert(this->current_b==0);
	assert(ibands >= 0 && ibands < this->nbands);	
    assert(ibasis >= 0 && ibasis < this->nbasis);
    if (this->device != CpuDevice) {
        #if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::WARNING_QUIT("Psi operator ()", "GPU psi cannot fetch value by an overloaded operator ()!");
        #endif 
    }
	return this->psi_current[ibands * this->nbasis + ibasis];
}

template<typename T, typename Device>
T& Psi<T, Device>::operator()(const int ibasis)
{	
    assert(ibasis >= 0 && ibasis < this->nbasis);
    if (this->device != CpuDevice) {
        #if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::WARNING_QUIT("Psi operator ()", "GPU psi cannot fetch value by an overloaded operator ()!");
        #endif 
    }
	return this->psi_current[this->current_b * this->nbasis + ibasis];
}

template<typename T, typename Device>
const T& Psi<T, Device>::operator()(const int ibasis) const
{	
    assert(ibasis >= 0 && ibasis < this->nbasis);
    if (this->device != CpuDevice) {
        #if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::WARNING_QUIT("Psi operator ()", "GPU psi cannot fetch value by an overloaded operator ()!");
        #endif 
    }
	return this->psi_current[this->current_b * this->nbasis + ibasis];
}

template<typename T, typename Device>
int Psi<T, Device>::get_current_k() const
{	
    return this->current_k;
}

template<typename T, typename Device>
int Psi<T, Device>::get_current_b() const
{	
    return this->current_b;
}

template<typename T, typename Device>
int Psi<T, Device>::get_current_nbas() const
{	
    return this->current_nbasis;
}


template<typename T, typename Device>
const int& Psi<T, Device>::get_ngk(const int ik_in) const
{	
    return this->ngk[ik_in];
}

template<typename T, typename Device>
void Psi<T, Device>::zero_out()
{	
    // this->psi.assign(this->psi.size(), T(0));
    set_memory_op()(this->ctx, this->psi, 0, this->size());
}


template<typename T, typename Device>
std::tuple<const T*, int> Psi<T, Device>::to_range(const Range& range) const
{
    int index_1_in = range.index_1;
    //mem_saver=1 case, only k==0 memory space is avaliable
    if (index_1_in >0 & this->nk == 1) {
        index_1_in = 0;
    }
    if (range.k_first != this->k_first 
        || index_1_in < 0 
        || range.range_1 < 0 
        || range.range_2<range.range_1
        || (range.k_first && range.range_2 >= this->nbands)
        || (!range.k_first && (range.range_2 >= this->nk || range.index_1>=this->nbands))) 
    {
        return std::tuple<const T*, int>(nullptr, 0);
    }
    else {
        const T* p = &this->psi[(index_1_in * this->nbands + range.range_1) * this->nbasis];
        int m = (range.range_2 - range.range_1 + 1) * this->npol;
        return std::tuple<const T*, int>(p, m);
    }
}

// only iterative diagonaliztion need initialization of Psi
void initialize(Psi<std::complex<double>> &psi)
{
    return;
}

void initialize(Psi<double> &psi)
{
    return;
}


namespace psi {
template class Psi<std::complex<double>, DEVICE_CPU>;
template class Psi<double, DEVICE_CPU>;
template Psi<std::complex<double>, DEVICE_CPU>::Psi<std::complex<double>, DEVICE_CPU>(const Psi<std::complex<double>, DEVICE_CPU>&);
template Psi<double, DEVICE_CPU>::Psi<double, DEVICE_CPU>(const Psi<double, DEVICE_CPU>&);
#if ((defined __CUDA) || (defined __ROCM))
template class Psi<std::complex<double>, DEVICE_GPU>;
template class Psi<double, DEVICE_GPU>;
template Psi<std::complex<double>, DEVICE_GPU>::Psi<std::complex<double>, DEVICE_GPU>(const Psi<std::complex<double>, DEVICE_GPU>&);
template Psi<double, DEVICE_GPU>::Psi<double, DEVICE_GPU>(const Psi<double, DEVICE_GPU>&);
template Psi<std::complex<double>, DEVICE_GPU>::Psi<std::complex<double>, DEVICE_CPU>(const Psi<std::complex<double>, DEVICE_CPU>&);
template Psi<double, DEVICE_GPU>::Psi<double, DEVICE_CPU>(const Psi<double, DEVICE_CPU>&);
template Psi<std::complex<double>, DEVICE_CPU>::Psi<std::complex<double>, DEVICE_GPU>(const Psi<std::complex<double>, DEVICE_GPU>&);
template Psi<double, DEVICE_CPU>::Psi<double, DEVICE_GPU>(const Psi<double, DEVICE_GPU>&);
#endif
} // namespace psi