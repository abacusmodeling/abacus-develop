#include "NVE.h"
#include "MD_func.h"

NVE::NVE(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{}

NVE::~NVE(){}

void NVE::setup()
{
    ModuleBase::TITLE("NVE", "setup");
    ModuleBase::timer::tick("NVE", "setup");

    MD_func::force_virial(mdp, ucell, potential, force, stress);
    double kenetic=MD_func::GetAtomKE(ucell.nat, vel, allmass);
    double tempNow = 2*kenetic/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU;
    MD_func::outStress(ucell, stress, kenetic);




    ModuleBase::timer::tick("NVE", "setup");
}

void NVE::first_half()
{
    ModuleBase::TITLE("NVE", "first_half");
    ModuleBase::timer::tick("NVE", "first_half");

    //MD_func::force_stress(mdp, ucell, potential, force, stress);

    ModuleBase::timer::tick("NVE", "first_half");
}

void NVE::second_half()
{
    ModuleBase::TITLE("NVE", "second_half");
    ModuleBase::timer::tick("NVE", "second_half");

    //MD_func::force_stress(mdp, ucell, potential, force, stress);

    ModuleBase::timer::tick("NVE", "second_half");
}
