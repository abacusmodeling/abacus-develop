#include "./esolver_ks_pw.h"
#include "../src_pw/sto_wf.h"
#include "../src_pw/sto_iter.h"
#include "../src_pw/sto_hchi.h"

namespace ModuleESolver
{

class ESolver_SDFT_PW: public ESolver_KS_PW
{
public:
    ESolver_SDFT_PW();
    ~ESolver_SDFT_PW();
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix& force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;
public:
    Stochastic_WF stowf;
    Stochastic_Iter stoiter;

protected:
    virtual void beforescf(const int istep) override; 
    // virtual void eachiterinit(int iter) override; 
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void eachiterfinish(const int iter) override; 
    virtual void afterscf() override;
    virtual void postprocess() override;

public:
    //Stochastic Kubo-Greenwood
    void sKG(const int nche_KG, const double fwhmin, const double wcut, 
             const double dw_in, const int times);
    void sKG_new(const int nche_KG, const double fwhmin, const double wcut, 
             const double dw_in, const int times);

private:
    int nche_sto; //norder of Chebyshev


};

}

