#include "surchem.h"

namespace GlobalC
{
  surchem solvent_model;
}

surchem::surchem()
{
    TOTN_real = new double[1];
    delta_phi = new double[1];
    Vcav = ModuleBase::matrix();
    Vel = ModuleBase::matrix();
    qs = 0;
}

void surchem::allocate(const int &nrxx, const int &nspin)
{
    assert(nrxx > 0);
    assert(nspin > 0);

    delete[] TOTN_real;
    delete[] delta_phi;

    TOTN_real = new double[nrxx];
    delta_phi = new double[nrxx];
    Vcav.create(nspin, nrxx);
    Vel.create(nspin, nrxx);

    ModuleBase::GlobalFunc::ZEROS(delta_phi, nrxx);
    ModuleBase::GlobalFunc::ZEROS(TOTN_real, nrxx);
    return;
}

surchem::~surchem()
{
    delete[] TOTN_real;
    delete[] delta_phi;
}