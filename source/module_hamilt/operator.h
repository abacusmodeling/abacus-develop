#ifndef __OPERATOR
#define __OPERATOR

#include<complex>
#include "module_psi/psi.h"
#include "module_base/global_function.h"

namespace hamilt
{

class Operator
{
    public:
    Operator(){};
    ~Operator()
    { 
        if(this->hpsi != nullptr) delete this->hpsi;
        this->release();
        Operator* last = this->next_op;
        if(last != nullptr) 
        {
            last->release();
            Operator* node_delete = last;
            last = last->next_op;
            node_delete->next_op = nullptr;
            delete node_delete;
        } 
    }
    virtual void release(){return;}

    typedef std::tuple<const psi::Psi<std::complex<double>>*, const psi::Range> hpsi_info;
    virtual hpsi_info hPsi(const hpsi_info& input)const {return hpsi_info(nullptr, 0);}

    virtual void act(std::complex<double> *hk_matrix)const {return;}
    virtual void act(double *hk_matrix)const {return;}

    virtual void init(const int ik_in)
    {
        this->ik = ik_in;
        if(this->next_op != nullptr)
        {
            this->next_op->init(ik_in);
        }
    }

    virtual void add(Operator* next)
    {
        if(next==nullptr) return;
        if(next->next_op != nullptr) this->add(next->next_op);
        Operator* last = this;
        while(last->next_op != nullptr)
        {
            if(next->cal_type==last->cal_type)
            {
                last->add(next);
                return;
            }
            last = last->next_op;
        }
        last->next_op = next; 
    }

    protected:
    int ik = 0;

    //calculation type, only different type can be in main chain table 
    int cal_type = 0;
    Operator* next_op = nullptr;

    //if this Operator is first node in chain table, hpsi would not be empty
    mutable psi::Psi<std::complex<double>>* hpsi = nullptr;

    std::complex<double>* get_hpsi(const hpsi_info& info)const
    {
        const int nbands_range = (std::get<1>(info).range_2 - std::get<1>(info).range_1 + 1);
        if(this->hpsi == nullptr)
        {
            delete this->hpsi;
        }
        this->hpsi = new psi::Psi<std::complex<double>>(std::get<0>(info)[0], 1, nbands_range);
        
        std::complex<double>* pointer_hpsi = this->hpsi->get_pointer();
        size_t total_hpsi_size = nbands_range * this->hpsi->get_nbasis();
        ModuleBase::GlobalFunc::ZEROS(pointer_hpsi, total_hpsi_size);
        return pointer_hpsi;
    }
};

}//end namespace hamilt

#endif