#include "esolver_fp.h"
#include "../module_base/global_variable.h"
namespace ModuleESolver
{
    void ESolver_FP::Init(Input& inp, UnitCell_pseudo& cell)
    {
        this->pw_rho.initgrids(cell.lat0, cell.latvec, 4 * inp.ecutwfc, GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
        this->pw_rho.initparameters(false, 4 * inp.ecutwfc);
        this->pw_rho.setuptransform();
        this->pw_rho.collect_local_pw(); 
    }
}