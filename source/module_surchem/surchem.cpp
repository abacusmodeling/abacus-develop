#include "surchem.h"

namespace GlobalC
{
  surchem solvent_model;
}

surchem::surchem()
{
    TOTN_real = new double[1];
    delta_phi = new double[1];
}

void surchem::allocate(const int &nrxx)
{
    delete[] TOTN_real;
    delete[] delta_phi;
    TOTN_real = new double[nrxx];
    delta_phi = new double[nrxx];
}

surchem::~surchem()
{
    delete[] TOTN_real;
    delete[] delta_phi;
}