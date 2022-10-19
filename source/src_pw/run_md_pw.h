//=======================
#ifndef RUN_MD_PW_H
#define RUN_MD_PW_H

#include "../src_pw/charge_extra.h"
#include "module_esolver/esolver.h"

class Run_MD_PW
{
public:
    Run_MD_PW();
    ~Run_MD_PW();

    void md_ions_pw(ModuleESolver::ESolver *p_esolver);

private:
    Charge_Extra CE;
    bool cellchange;
};

#endif