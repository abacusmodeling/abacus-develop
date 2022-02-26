#include "Cell_PW.h"
#include "ions.h"
#include "../module_base/timer.h"

void Cell_PW::opt_cells_pw(ModuleEnSover::En_Solver *p_ensolver)
{
    ModuleBase::TITLE("Cell_PW", "opt_cells_pw");
    ModuleBase::timer::tick("Cell_PW", "opt_cells_pw");


    // ion optimization begins
    // electron density optimization is included in ion optimization

    Ions ions;
    ions.opt_ions_pw(p_ensolver);
    
    ModuleBase::timer::tick("Cell_PW", "opt_cells_pw");
}
