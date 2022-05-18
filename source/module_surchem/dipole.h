#ifndef DIPOLE_H
#define DIPOLE_H

#include "../src_pw/pw_basis.h"
#include "../module_cell/unitcell.h"

class Dipole
{
public:
    Dipole();
    ~Dipole();

    static ModuleBase::matrix add_efield(const UnitCell &cell, 
                                            PW_Basis &pwb, 
                                            const int &nspin, 
                                            const double *const *const rho);

    static double cal_elec_dipole(const UnitCell &cell, 
                                PW_Basis &pwb, 
                                const int &nspin, 
                                const double *const *const rho,
                                const double &h_inv);

    static double cal_ion_dipole(const UnitCell &cell, const double &h_inv);

    static double saw_function(const double &a, const double &b, const double &x);

    static void compute_force(const UnitCell &cell, ModuleBase::matrix &fdip);



    static double etotefield;           // dipole energy
    static double tot_dipole;           // total dipole
    static int edir;                    // 0, 1, 2 denotes x, y, z direction for dipole correction
    static double emaxpos;              // the maximum position of the saw function
    static double eopreg;               // the decrease region length of the saw function
    static double eamp;                 // field amplitude (in a.u.) (1 a.u. = 51.44 10^10 V/m)
    static bool first;                  // efield only needs to be added on the first iteration, if dipfield is not used.
    static double bvec[3];
    static double bmod;
};





#endif