#include "surchem.h"

namespace GlobalC
{
surchem solvent_model;
}

surchem::surchem()
{
    TOTN_real = nullptr;
    delta_phi = nullptr;
    epspot = nullptr;
    comp_real = nullptr;
    phi_comp_R = nullptr;
    Vcav = ModuleBase::matrix();
    Vel = ModuleBase::matrix();
    qs = 0;
}

void surchem::allocate(const int &nrxx, const int &nspin)
{
    assert(nrxx >= 0);
    assert(nspin > 0);

    delete[] TOTN_real;
    delete[] delta_phi;
    delete[] epspot;
    delete[] comp_real;
    delete[] phi_comp_R;
    if (nrxx > 0)
    {
        TOTN_real = new double[nrxx];
        delta_phi = new double[nrxx];
        epspot = new double[nrxx];
        comp_real = new double[nrxx];
        phi_comp_R = new double[nrxx];
    }
    else
    {
        TOTN_real = nullptr;
        delta_phi = nullptr;
        epspot = nullptr;
        comp_real = nullptr;
        phi_comp_R = nullptr;
    }
    Vcav.create(nspin, nrxx);
    Vel.create(nspin, nrxx);

    ModuleBase::GlobalFunc::ZEROS(delta_phi, nrxx);
    ModuleBase::GlobalFunc::ZEROS(TOTN_real, nrxx);
    ModuleBase::GlobalFunc::ZEROS(epspot, nrxx);
    ModuleBase::GlobalFunc::ZEROS(comp_real, nrxx);
    ModuleBase::GlobalFunc::ZEROS(phi_comp_R, nrxx);
    return;
}

surchem::~surchem()
{
    delete[] TOTN_real;
    delete[] delta_phi;
    delete[] epspot;
    delete[] comp_real;
    delete[] phi_comp_R;
}