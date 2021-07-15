#ifndef RUN_MD_CLASSIC_H
#define RUN_MD_CLASSIC_H

#include "../src_pw/tools.h"
#include "../module_neighbor/sltk_grid_driver.h"
#include "../module_cell/unitcell_pseudo.h"

class Run_MD_CLASSIC
{
public:
    Run_MD_CLASSIC();
    ~Run_MD_CLASSIC();

    Grid_Driver grid_neigh;
    UnitCell_pseudo ucell_c;
    

    void classic_md_line(void);
    void md_cells_classic(void);
    void md_allocate_ions(void);
    void update_pos_classic(void);


private:
    int istep;
    double* pos_old1;
	double* pos_old2;
	double* pos_now;
	double* pos_next;
    int pos_dim;
   
};

#endif