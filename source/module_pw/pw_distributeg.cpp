#include "pw_basis.h"
#include "../module_base/tool_quit.h"

namespace ModulePW
{
//
//distribute plane waves to different cores
//Known: G, GT, GGT, nx, ny, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2ixy[is], ixy2istot[ixy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nstnz_per[ip], gg[ig], gcar[ig], gdirect[ig], nst, nstot
//
void PW_Basis::distribute_g()
{
    if(this->distribution_type == 1)
    {
        this->distribution_method1();
    }
    else if(this->distribution_type == 2)
    {
        this->distribution_method2();
    }
    else
    {
        ModuleBase::WARNING_QUIT("divide", "No such division type.");
    }
    return;
}
}
