#include "module_hamilt_general/operator.h"

using namespace hamilt;


template<typename FPTYPE, typename Device>
Operator<FPTYPE, Device>::Operator(){}

template<typename FPTYPE, typename Device>
Operator<FPTYPE, Device>::~Operator() 
{
    if(this->hpsi != nullptr) delete this->hpsi;
    Operator* last = this->next_op;
    Operator* last_sub = this->next_sub_op;
    while(last != nullptr || last_sub != nullptr)
    {
        if(last_sub != nullptr)
        {//delete sub_chain first
            Operator* node_delete = last_sub;
            last_sub = last_sub->next_sub_op;
            node_delete->next_sub_op = nullptr;
            delete node_delete;
        }
        else
        {//delete main chain if sub_chain is deleted
            Operator* node_delete = last;
            last_sub = last->next_sub_op;
            node_delete->next_sub_op = nullptr;
            last = last->next_op;
            node_delete->next_op = nullptr;
            delete node_delete;
        }
    }
}

template<typename FPTYPE, typename Device>
typename Operator<FPTYPE, Device>::hpsi_info Operator<FPTYPE, Device>::hPsi(hpsi_info&) const 
{
    ModuleBase::WARNING_QUIT("Operator::hPsi", "hPsi error!");
    return hpsi_info(nullptr, 0, nullptr);
}

template<typename FPTYPE, typename Device>
void Operator<FPTYPE, Device>::init(const int ik_in) 
{
    this->ik = ik_in;
    if(this->next_op != nullptr) {
        this->next_op->init(ik_in);
    }
}

template<typename FPTYPE, typename Device>
void Operator<FPTYPE, Device>::add(Operator* next) 
{
    if(next==nullptr) return;
    next->is_first_node = false;
    if(next->next_op != nullptr) this->add(next->next_op);
    Operator* last = this;
    //loop to end of the chain
    while(last->next_op != nullptr)
    {
        if(next->cal_type==last->cal_type)
        {
            break;
        }
        last = last->next_op;
    }
    if(next->cal_type == last->cal_type)
    {
        //insert next to sub chain of current node
        Operator* sub_last = last;
        while(sub_last->next_sub_op != nullptr)
        {
            sub_last = sub_last->next_sub_op;
        }
        sub_last->next_sub_op = next;
        return;
    }
    else
    {
        last->next_op = next;
    }
}

template<typename FPTYPE, typename Device>
FPTYPE* Operator<FPTYPE, Device>::get_hpsi(const hpsi_info& info) const
{
    const int nbands_range = (std::get<1>(info).range_2 - std::get<1>(info).range_1 + 1);
    //in_place call of hPsi, hpsi inputs as new psi, 
    //create a new hpsi and delete old hpsi later
    FPTYPE* hpsi_pointer = std::get<2>(info);
    const FPTYPE* psi_pointer = std::get<0>(info)->get_pointer();
    if(this->hpsi != nullptr) 
    {
        delete this->hpsi;
        this->hpsi = nullptr;
    }
    if(!hpsi_pointer)
    {
        ModuleBase::WARNING_QUIT("Operator::hPsi", "hpsi_pointer can not be nullptr");
    }
    else if(hpsi_pointer == psi_pointer)
    {
        this->in_place = true;
        this->hpsi = new psi::Psi<FPTYPE, Device>(std::get<0>(info)[0], 1, nbands_range);
    }
    else
    {
        this->in_place = false;
        this->hpsi = new psi::Psi<FPTYPE, Device>(hpsi_pointer, std::get<0>(info)[0], 1, nbands_range);
    }
    
    hpsi_pointer = this->hpsi->get_pointer();
    size_t total_hpsi_size = nbands_range * this->hpsi->get_nbasis();
    // ModuleBase::GlobalFunc::ZEROS(hpsi_pointer, total_hpsi_size);
    // denghui replaced at 20221104
    set_memory_op()(this->ctx, hpsi_pointer, 0, total_hpsi_size);
    return hpsi_pointer;
}


namespace hamilt {
template class Operator<float, psi::DEVICE_CPU>;
template class Operator<std::complex<float>, psi::DEVICE_CPU>;
template class Operator<double, psi::DEVICE_CPU>;
template class Operator<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Operator<float, psi::DEVICE_GPU>;
template class Operator<std::complex<float>, psi::DEVICE_GPU>;
template class Operator<double, psi::DEVICE_GPU>;
template class Operator<std::complex<double>, psi::DEVICE_GPU>;
#endif
}
