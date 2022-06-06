#include "FORCE_gamma.h"
#include "../src_pw/global.h"
#include "../module_base/timer.h"

void Force_LCAO_gamma::cal_fvl_dphi(
	double*** DM_in,
	const bool isforce, 
    const bool isstress,
    ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& svl_dphi)
{   
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
    int istep = 1;
    GlobalC::pot.init_pot(istep, GlobalC::sf.strucFac);
    fvl_dphi.zero_out();
    svl_dphi.zero_out();
    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;
        for(int ir=0; ir<GlobalC::rhopw->nrxx; ++ir)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        Gint_inout inout(DM_in, GlobalC::pot.vr_eff1, isforce, isstress, &fvl_dphi, &svl_dphi, Gint_Tools::job_type::force);

        this->UHM->GG.cal_gint(&inout);
    }

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= -GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
}