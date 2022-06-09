#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "./esolver_ks_lcao.h"

#include "src_lcao/record_adj.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_hamilt.h"
#include "module_orbital/ORB_control.h"

namespace ModuleESolver
{

    class ESolver_KS_LCAO_TDDFT : ESolver_KS_LCAO
    {
    public:
        ESolver_KS_LCAO_TDDFT();
        ~ESolver_KS_LCAO_TDDFT();
        void Init(Input& inp, UnitCell_pseudo& cell) override;

    protected:
        virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
        virtual void eachiterinit(const int istep, const int iter) override;
        virtual void updatepot(const int istep, const int iter) override;
        void cal_edm_tddft();
    };



}
#endif
