#ifndef DIPOLE_H
#define DIPOLE_H

#include "pw_basis.h"
#include "../module_cell/unitcell.h"

class Dipole
{
public:
    Dipole();
    ~Dipole();

    static double dipole_density(const UnitCell &cell, PW_Basis &pwb, const int &nspin, const double *const *const rho);
    static ModuleBase::matrix v_dipole(const UnitCell &cell, PW_Basis &pwb, const int &nspin, const double *const *const rho);



    static double dipole_energy;        // dipole energy
    static int dir;                     // 0, 1, 2 denotes x, y, z direction for dipole correction
};





#endif