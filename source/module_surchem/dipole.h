#ifndef DIPOLE_H
#define DIPOLE_H

#include "../src_pw/pw_basis.h"
#include "../module_cell/unitcell.h"

class Dipole
{
public:
    Dipole();
    ~Dipole();

    static double dipole(const UnitCell &cell, 
                            PW_Basis &pwb, 
                            const double *const *const rho, 
                            complex<double> *TOTN);

    static ModuleBase::matrix v_dipole(const UnitCell &cell, 
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



    static double dipole_energy;        // dipole energy
    static int dir;                     // 0, 1, 2 denotes x, y, z direction for dipole correction
    static double max_pos;              // the maximum position of the saw function
    static double de_reg;               // the decrease region length of the saw function
};





#endif