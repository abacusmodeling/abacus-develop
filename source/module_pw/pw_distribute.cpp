#include "pw_basis.h"

//
//distribute plane waves to different processors
//Known: G, GT, GGT, nx, ny, nz, poolnproc, poolrank, ggecut
//output: ig2fft[ig], is2ir[is], gg[ig], gcar[ig], gdirect[ig], nst
//
void PW_Basis::distribute_g()
{
    if(this->divide_type == 1)
    {
        this->distribution_method1();
    }
    else
    {
        ModuleBase::WARNING_QUIT("divide", "No such division type.");
    }
    return;
}

void PW_Basis::distribution_method1()
{
    return;
}

