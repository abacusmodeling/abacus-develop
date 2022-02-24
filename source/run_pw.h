//==========================================================
// AUTHOR : mohan
// DATE : 2021-02-01
//==========================================================
#ifndef RUN_PW_H
#define RUN_PW_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "input.h"
#include "module_ensolver/en_solver.h"

class Run_pw
{

	public:

    Run_pw();
    ~Run_pw();

	// perform plane wave basis calculations
    static void plane_wave_line(ModuleEnSover::En_Solver *p_ensolver);

};

#endif
