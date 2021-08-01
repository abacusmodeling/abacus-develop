#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H

#include "../module_base/vector3.h"
#include "../module_base/matrix.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../module_neighbor/sltk_grid_driver.h"

class LJ_potential
{
public:

    LJ_potential();
    ~LJ_potential();

    double Lennard_Jones(UnitCell_pseudo &ucell_c, Grid_Driver &grid_neigh, Vector3<double> *force, matrix &stress);
    double LJ_energy(const double d);
    Vector3<double> LJ_force(const double d, const Vector3<double> dr);

};

#endif