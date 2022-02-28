#ifndef KS_SCF_LCAO_H
#define KS_SCF_LCAO_H
#include "../ks_scf.h"

#include "src_lcao/LOOP_elec.h"

namespace ModuleEnSover
{

class KS_SCF_LCAO: public KS_SCF
{
public:
    KS_SCF_LCAO()
    {
        tag = "KS_SCF_LCAO";
    }
    void Init(Input& inp, UnitCell_pseudo& cell) override;
    
    void Run(int istep,
        Record_adj& ra,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& lowf,
        LCAO_Hamilt& uhm) override;
    void Run(int istep, UnitCell_pseudo& cell) override {};
    
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;

private:
    LOOP_elec LOE;
};

///Basis_lcao
///ORB_control orb_con;


}
#endif