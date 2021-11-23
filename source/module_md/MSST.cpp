#include "MSST.h"
#include "MD_func.h"

MSST::MSST(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    std::cout << "MSST" << std::endl;

    mdp.Qmass = mdp.Qmass / pow(ModuleBase::ANGSTROM_AU, 4) / pow(ModuleBase::AU_to_MASS, 2);
    mdp.velocity = mdp.velocity * ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;
    mdp.viscosity = mdp.viscosity / ModuleBase::AU_to_MASS / ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;

    

}

MSST::~MSST(){}

void MSST::setup(){}

void MSST::first_half(){}

void MSST::second_half(){}

void MSST::outputMD(){}