#include "./esolver_ks_pw.h"
namespace ModuleESolver
{

class ESolver_SDFT_PW: public ESolver_KS_PW
{
public:
    ESolver_SDFT_PW();
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix &stress) override;

protected:
    virtual void beforescf() override; 
    // virtual void eachiterinit(int iter) override; 
    virtual void hamilt2density(int istep, int iter, double ethr) override;
    virtual void eachiterfinish(int iter, bool conv) override; 
    virtual void afterscf(bool) override;
};

}

