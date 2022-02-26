#ifndef KS_SCF_PW_H
#define KS_SCF_PW_H
#include "../ks_scf.h"
// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"

#include "../../../../src_pw/electrons.h"

namespace ModuleEnSover
{

class KS_SCF_PW: public KS_SCF
{
public:
    KS_SCF_PW()
    {
        tag = "KS_SCF_PW";
    }
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void Run(int istep, UnitCell_pseudo &cell) override;
    void cal_Energy(energy &en) override; 
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix &stress) override;

//--------------temporary----------------------------
    Electrons elec;
//---------------------------------------------------

    // KS_SCF_PW(bool use_sdft)
    // {
    //     if(use_sdft)    
    //     {
    //         this->ph2e=new H2E_SDFT();
    //     }
    //     else           
    //     {
    //         this->ph2e=new H2E_PW();
    //     }
    //     this->pes= new Estate_PW();
    //     this->phamilt=new Hamilt_PW();
    // }
    
    // Basis_PW basis_pw;
    // Init(Inputs &inp, Cell &cel)
    // {
        
    //     basis_pw.init(inp, cel);

    //     pes->init(inp, cel, basis_pw); 
        
    //     phamilt->init(bas); 
    //     phamilt->initpot(cel, pes); 
        
    //     ph2e->init(h, pes); 
    // }
};
}
#endif