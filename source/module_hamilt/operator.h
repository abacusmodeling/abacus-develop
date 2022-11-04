#ifndef OPERATOR_H
#define OPERATOR_H

#include<complex>
#include "module_psi/psi.h"
#include "module_base/global_function.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

enum calculation_type
{
    no,
    pw_ekinetic, 
    pw_nonlocal,
    pw_veff,
    pw_meta,
    lcao_fixed,
    lcao_gint,
    lcao_deepks,
    lcao_exx,
    lcao_dftu,
};

// Basic class for operator module, 
// it is designed for "O|psi>" and "<psi|O|psi>"
// Operator "O" might have several different types, which should be calculated one by one.
// In basic class , function add() is designed for combine all operators together with a chain. 
template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Operator
{
    public:
    Operator();
    virtual ~Operator();

    //this is the core function for Operator
    // do H|psi> from input |psi> , 

    // output of hpsi would be first member of the returned tuple 
    typedef std::tuple<const psi::Psi<FPTYPE, Device>*, const psi::Range, FPTYPE*> hpsi_info;
    virtual hpsi_info hPsi(hpsi_info& input)const;

    virtual void init(const int ik_in);

    virtual void add(Operator* next);

    #if ((defined __CUDA) || (defined __ROCM))
    typedef std::tuple<const psi::Psi<FPTYPE, psi::DEVICE_GPU>*, const psi::Range, FPTYPE*> hpsi_info_gpu;
    virtual hpsi_info_gpu hPsi_gpu(hpsi_info_gpu& input) const; 
    #endif // ((defined __CUDA) || (defined __ROCM))


    protected:
    int ik = 0;

    mutable bool in_place = false;

    //calculation type, only different type can be in main chain table 
    enum calculation_type cal_type;
    Operator* next_op = nullptr;
    Operator* next_sub_op = nullptr;
    bool is_first_node = true;

    //if this Operator is first node in chain table, hpsi would not be empty
    mutable psi::Psi<FPTYPE, Device>* hpsi = nullptr;

    /*This function would analyze hpsi_info and choose how to arrange hpsi storage
    In hpsi_info, if the third parameter hpsi_pointer is set, which indicates memory of hpsi is arranged by developer;
    if hpsi_pointer is not set(nullptr), which indicates memory of hpsi is arranged by Operator, this case is rare. 
    two cases would occurred:
    1. hpsi_pointer != nullptr && psi_pointer == hpsi_pointer , psi would be replaced by hpsi, hpsi need a temporary memory
    2. hpsi_pointer != nullptr && psi_pointer != hpsi_pointer , this is the commonly case 
    */
    FPTYPE* get_hpsi(const hpsi_info& info)const;

};

}//end namespace hamilt

#endif