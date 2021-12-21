#include "Langevin.h"
#include "MD_func.h"

Langevin::Langevin(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    // convert to a.u. unit
    mdp.damp /= ModuleBase::AU_to_FS;
}

Langevin::~Langevin(){}

void Langevin::setup()
{
    ModuleBase::TITLE("Langevin", "setup");
    ModuleBase::timer::tick("Langevin", "setup");
    
    Verlet::setup();

    post_force();

    ModuleBase::timer::tick("Langevin", "setup");
}

void Langevin::first_half()
{
    ModuleBase::TITLE("Langevin", "first_half");
    ModuleBase::timer::tick("Langevin", "first_half");

    Verlet::first_half();

    ModuleBase::timer::tick("Langevin", "first_half");
}

void Langevin::second_half()
{
    ModuleBase::TITLE("Langevin", "second_half");
    ModuleBase::timer::tick("Langevin", "second_half");

    post_force();

    Verlet::second_half();

    ModuleBase::timer::tick("Langevin", "second_half");
}

void Langevin::outputMD()
{
    Verlet::outputMD();
}

void Langevin::write_restart()
{
    Verlet::write_restart();
}

void Langevin::restart()
{
    Verlet::restart();
}

void Langevin::post_force()
{
    temp_target();

    for(int i=0; i<ucell.nat; ++i)
    {
        force[i] -= allmass[i] * vel[i] / mdp.damp;
        for(int j=0; j<3; ++j)
        {
            force[i][j] += sqrt(24.0 * t_target * allmass[i] / mdp.damp / mdp.dt) * (rand()/double(RAND_MAX) - 0.5);
        }
    }
}

void Langevin::temp_target()
{
    double delta = (double)(step_ + step_rst_) / GlobalV::NSTEP;
    t_target = mdp.tfirst + delta * (mdp.tlast - mdp.tfirst);
}