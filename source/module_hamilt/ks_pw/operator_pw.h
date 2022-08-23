#ifndef __OPERATORPW
#define __OPERATORPW
#include"module_hamilt/operator.h"
#include "module_base/timer.h"

namespace hamilt
{

class OperatorPW : public Operator<std::complex<double>>
{
    public:
    virtual ~OperatorPW(){};
    
    //in PW code, different operators donate hPsi independently
    //run this->act function for the first operator and run all act() for other nodes in chain table  
    virtual hpsi_info hPsi(const hpsi_info& input)const
    {
        ModuleBase::timer::tick("OperatorPW", "hPsi");
        std::tuple<const std::complex<double>*, int> psi_info = std::get<0>(input)->to_range(std::get<1>(input));
        int n_npwx = std::get<1>(psi_info); 

        std::complex<double> *tmhpsi = this->get_hpsi(input);
        const std::complex<double> *tmpsi_in = std::get<0>(psi_info);
        //if range in hpsi_info is illegal, the first return of to_range() would be nullptr
        if(tmpsi_in == nullptr)
        {
            ModuleBase::WARNING_QUIT("OperatorPW", "please choose correct range of psi for hPsi()!");
        }

        this->act(std::get<0>(input), n_npwx, tmpsi_in, tmhpsi);
        OperatorPW* node((OperatorPW*)this->next_op);
        while(node != nullptr)
        {
            node->act(std::get<0>(input), n_npwx, tmpsi_in, tmhpsi);
            node = (OperatorPW*)(node->next_op);
        }

        //during recursive call of hPsi, delete the input psi
        if(this->recursive) delete std::get<0>(input);

        ModuleBase::timer::tick("OperatorPW", "hPsi");
        
        return hpsi_info(this->hpsi, psi::Range(1, 0, 0, n_npwx/std::get<0>(input)->npol));
    }
    
    //main function which should be modified in Operator for PW base
    //do operation : |hpsi_choosed> = V|psi_choosed>
    //V is the target operator act on choosed psi, the consequence should be added to choosed hpsi
    virtual void act
    (
        const psi::Psi<std::complex<double>> *psi_in, 
        const int n_npwx, 
        const std::complex<double>* tmpsi_in, 
        std::complex<double>* tmhpsi
    )const 
    {
        return;
    }


};

}//end namespace hamilt

#endif