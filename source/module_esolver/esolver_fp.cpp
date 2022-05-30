#include "esolver_fp.h"
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
namespace ModuleESolver
{   ESolver_FP::ESolver_FP()
    {
        // pw_rho = new ModuleBase::PW_Basis();
        //temporary, it will be removed
        pw_rho = new ModulePW::PW_Basis_Big(); 
        ModulePW::PW_Basis_Big* tmp = static_cast<ModulePW::PW_Basis_Big*>(pw_rho);
        tmp->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    }
    ESolver_FP::~ESolver_FP()
    {
        if(pw_rho!=NULL)    delete pw_rho;
    }
    void ESolver_FP::Init(Input& inp, UnitCell_pseudo& cell)
    {
        this->pw_rho->initgrids(cell.lat0, cell.latvec, inp.ecutrho, GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
        this->pw_rho->initparameters(false, inp.ecutrho);
        this->pw_rho->setuptransform();
        this->pw_rho->collect_local_pw(); 
        this->pw_rho->collect_uniqgg();
        GlobalC::rhopw = this->pw_rho; //Temporary
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of |g|", this->pw_rho->ngg);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max |g|", this->pw_rho->gg_uniq[ this->pw_rho->ngg-1]);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"min |g|", this->pw_rho->gg_uniq[0]);

    }
}