#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

void Gint_Gamma::cal_force(double*** DM_in, const double*const vlocal, 
        ModuleBase::matrix& force, ModuleBase::matrix& stress, 
        const bool is_force, const bool is_stress)
{
    ModuleBase::TITLE("Grid_Integral","cal_force");
    ModuleBase::timer::tick("Gint_Gamma","cal_force");
    if(!is_force && !is_stress)
    {
        ModuleBase::timer::tick("Gint_Gamma","cal_force");
        return;
    }
    this->max_size = GlobalC::GridT.max_atom;

	if(max_size)
    {
        const int nbx = GlobalC::GridT.nbx;
        const int nby = GlobalC::GridT.nby;
        const int nbz_start = GlobalC::GridT.nbzp_start;
        const int nbz = GlobalC::GridT.nbzp;
        const double dv = GlobalC::ucell.omega/GlobalC::rhopw->nxyz;
    
        const int ncyz = GlobalC::rhopw->ny*GlobalC::rhopw->nplane; // mohan add 2012-03-25

        for (int i=0; i<nbx; i++)
        {
            const int ibx = i*GlobalC::bigpw->bx;
            for (int j=0; j<nby; j++)
            {
                const int jby = j*GlobalC::bigpw->by;
                for (int k=nbz_start; k<nbz_start+nbz; k++)
                {
                    const int kbz = k*GlobalC::bigpw->bz-GlobalC::rhopw->startz_current;
    
                    const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
    
                    // get the value: how many atoms has orbital value on this grid.
                    const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
                    if(na_grid==0) continue;

                    // it's a uniform grid to save orbital values, so the delta_r is a constant.
                    const double delta_r = GlobalC::ORB.dr_uniform;
                    
                    int *block_iw = new int[na_grid];
                    int *block_index = new int[na_grid+1];
                    int *block_size = new int[na_grid];
                    Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);

                    //------------------------------------------------------
                    // whether the atom-grid distance is larger than cutoff
                    //------------------------------------------------------
                    bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

                    // set up band matrix psir_ylm and psir_DM
                    const int LD_pool = max_size*GlobalC::ucell.nwmax;

                    Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::bigpw->bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_x(GlobalC::bigpw->bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_y(GlobalC::bigpw->bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_z(GlobalC::bigpw->bxyz, LD_pool);

                    Gint_Tools::cal_dpsir_ylm(
                        na_grid, grid_index, delta_r,
                        block_index, block_size, 
                        cal_flag,
                        psir_ylm.ptr_2D,
                        dpsir_ylm_x.ptr_2D,
                        dpsir_ylm_y.ptr_2D,
                        dpsir_ylm_z.ptr_2D
                    );

                    double* vldr3 = new double[GlobalC::bigpw->bxyz];
                    Gint_Tools::get_vldr3(vlocal, ncyz, ibx, jby, kbz, dv, vldr3);
                    const Gint_Tools::Array_Pool<double> psir_vlbr3    = Gint_Tools::get_psir_vlbr3(na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);
                    const Gint_Tools::Array_Pool<double> psir_vlbr3_DM = Gint_Tools::get_psir_vlbr3_DM(na_grid, LD_pool, block_iw, block_size, block_index, cal_flag, psir_vlbr3.ptr_2D, DM_in[GlobalV::CURRENT_SPIN]);

                    if(is_force)
                    {
                        this-> cal_meshball_force(grid_index, na_grid, 
                            block_size, block_index,
                            psir_vlbr3_DM.ptr_2D, 
                            dpsir_ylm_x.ptr_2D, 
                            dpsir_ylm_y.ptr_2D, 
                            dpsir_ylm_z.ptr_2D, 
                            force);
                    }
                    if(is_stress)
                    {
                        Gint_Tools::Array_Pool<double> dpsir_ylm_xx(GlobalC::bigpw->bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_xy(GlobalC::bigpw->bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_xz(GlobalC::bigpw->bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_yy(GlobalC::bigpw->bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_yz(GlobalC::bigpw->bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_zz(GlobalC::bigpw->bxyz, LD_pool);
                        Gint_Tools::cal_dpsirr_ylm(
                            na_grid, grid_index,
                            block_index, block_size, 
                            cal_flag,
                            dpsir_ylm_x.ptr_2D,
                            dpsir_ylm_y.ptr_2D,
                            dpsir_ylm_z.ptr_2D,
                            dpsir_ylm_xx.ptr_2D,
                            dpsir_ylm_xy.ptr_2D,
                            dpsir_ylm_xz.ptr_2D,
                            dpsir_ylm_yy.ptr_2D,
                            dpsir_ylm_yz.ptr_2D,
                            dpsir_ylm_zz.ptr_2D
                        );
                        this-> cal_meshball_stress(na_grid, block_index,
                            psir_vlbr3_DM.ptr_2D, 
                            dpsir_ylm_xx.ptr_2D, 
                            dpsir_ylm_xy.ptr_2D, 
                            dpsir_ylm_xz.ptr_2D,
                            dpsir_ylm_yy.ptr_2D, 
                            dpsir_ylm_yz.ptr_2D, 
                            dpsir_ylm_zz.ptr_2D,
                            stress);
                    }

                    delete[] vldr3;
                    delete[] block_iw;
                    delete[] block_index;
                    delete[] block_size;

                    for(int ib=0; ib<GlobalC::bigpw->bxyz; ++ib)
                    {
                        delete[] cal_flag[ib];
                    }
                    delete[] cal_flag;
                }//k
            }//j
        }//i
    }//max_size


    ModuleBase::timer::tick("Gint_Gamma","cal_force");
}

void Gint_Gamma::cal_meshball_force(
    const int grid_index,
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DM,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
    const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
    const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
    const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
    ModuleBase::matrix &force
)
{
    constexpr int inc=1;
    for(int ia1=0;ia1<na_grid;ia1++)
    {
        const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + ia1;
        const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom

        for(int ib=0;ib<GlobalC::bigpw->bxyz;ib++)
        {
            const double rx = ddot_(&block_size[ia1], &psir_vlbr3_DM[ib][block_index[ia1]], &inc, &dpsir_x[ib][block_index[ia1]], &inc);
            force(iat,0)+=rx*2.0;
            const double ry = ddot_(&block_size[ia1], &psir_vlbr3_DM[ib][block_index[ia1]], &inc, &dpsir_y[ib][block_index[ia1]], &inc);
            force(iat,1)+=ry*2.0;
            const double rz = ddot_(&block_size[ia1], &psir_vlbr3_DM[ib][block_index[ia1]], &inc, &dpsir_z[ib][block_index[ia1]], &inc);
            force(iat,2)+=rz*2.0;
          
        }
    }
    return;
}

void Gint_Gamma::cal_meshball_stress(
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DM,
    const double*const*const dpsir_xx,
    const double*const*const dpsir_xy,
    const double*const*const dpsir_xz,
    const double*const*const dpsir_yy,
    const double*const*const dpsir_yz,
    const double*const*const dpsir_zz,
    ModuleBase::matrix &stress
)
{
    constexpr int inc=1;
    for(int ib=0; ib<GlobalC::bigpw->bxyz; ++ib)
    {
        const double rxx = ddot_(&block_index[na_grid], psir_vlbr3_DM[ib], &inc, dpsir_xx[ib], &inc);
        stress(0,0)+=rxx*2.0;
        const double rxy = ddot_(&block_index[na_grid], psir_vlbr3_DM[ib], &inc, dpsir_xy[ib], &inc);
        stress(0,1)+=rxy*2.0;
        const double rxz = ddot_(&block_index[na_grid], psir_vlbr3_DM[ib], &inc, dpsir_xz[ib], &inc);
        stress(0,2)+=rxz*2.0;
        const double ryy = ddot_(&block_index[na_grid], psir_vlbr3_DM[ib], &inc, dpsir_yy[ib], &inc);
        stress(1,1)+=ryy*2.0;
        const double ryz = ddot_(&block_index[na_grid], psir_vlbr3_DM[ib], &inc, dpsir_yz[ib], &inc);
        stress(1,2)+=ryz*2.0;
        const double rzz = ddot_(&block_index[na_grid], psir_vlbr3_DM[ib], &inc, dpsir_zz[ib], &inc);
        stress(2,2)+=rzz*2.0;
    }
    return;
}