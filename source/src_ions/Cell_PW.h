//================================
// AUTHOR : zhengdy
// DATE : 2021-06-22
//================================
#ifndef CELL_PW_H
#define CELL_PW_H

#include "../src_pw/global.h"
#include "module_ensolver/en_solver.h"

class Cell_PW
{
public:
    Cell_PW(){};
    ~Cell_PW(){};
    
    void opt_cells_pw(ModuleEnSover::En_Solver *p_ensolver);
};

#endif