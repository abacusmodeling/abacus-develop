//=======================
#ifndef RUN_MD_PW_H
#define RUN_MD_PW_H

#include "../src_pw/tools.h"
#include "../src_pw/charge_extra.h"
#include "../src_ions/lattice_change_methods.h"
#include "../src_pw/global.h"

class Run_MD_PW
{
public:
    Run_MD_PW();
    ~Run_MD_PW();

    void md_ions_pw();
    void md_cells_pw();
    void md_force_virial(const int &istep,
        const int& numIon, 
        double &potential, 
        ModuleBase::Vector3<double>* force, 
        ModuleBase::matrix& virial);

private:
    Charge_Extra CE;
    Lattice_Change_Methods LCM;
};

#endif