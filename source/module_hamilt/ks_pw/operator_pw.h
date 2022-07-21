#ifndef __OPERATORPW
#define __OPERATORPW
#include"module_hamilt/operator.h"

namespace hamilt
{

class OperatorPW : public Operator
{
    public:
    
    //in PW code, different operators donate hPsi independently
    //run this->act function for the first operator and run all act() for other nodes in chain table  
    virtual hpsi_info hPsi(const hpsi_info& input)const
    {
        std::tuple<const std::complex<double>*, int> psi_info = std::get<0>(input)->to_range(std::get<1>(input));
        int n_npwx = std::get<1>(psi_info); 
        const int npwx = std::get<0>(input)->get_nbasis();

        std::complex<double> *tmhpsi = this->get_hpsi(input);
        const std::complex<double> *tmpsi_in = std::get<0>(psi_info);

        this->act(std::get<0>(input), n_npwx, tmpsi_in, tmhpsi);
        OperatorPW* node((OperatorPW*)this->next_op);
        while(node != nullptr)
        {
            node->act(std::get<0>(input), n_npwx, tmpsi_in, tmhpsi);
            node = (OperatorPW*)(node->next_op);
        }
        
        return hpsi_info(this->hpsi, std::get<1>(input));
    }
    
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