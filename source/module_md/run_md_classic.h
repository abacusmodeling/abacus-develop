#ifndef RUN_MD_CLASSIC_H
#define RUN_MD_CLASSIC_H

#include "../src_pw/tools.h"
#include "../module_cell/unitcell_pseudo.h"

class Run_MD_CLASSIC
{
public:
    Run_MD_CLASSIC();
    ~Run_MD_CLASSIC();

    UnitCell_pseudo ucell_c;
    

    void classic_md_line(void);
    void md_force_stress(double &potential);
    void md_allocate_ions(void);
    void update_pos_classic(void);

    Vector3<double> *force;  //force of each atom
	ModuleBase::matrix stress;           //stress for this lattice

private:
    int istep;
    double* pos_old1;
	double* pos_old2;
	double* pos_now;
	double* pos_next;
    int pos_dim;
};

#endif