#include "gint_k.h"
#include "../src_pw/global.h"
#include "global_fp.h" // mohan add 2021-01-30

#include "../module_base/ylm.h"
#include "../module_base/timer.h"

void Gint_k::cal_force_k(
	const bool isforce,
	const bool isstress,
	ModuleBase::matrix& fvl_dphi, 
	ModuleBase::matrix& svl_dphi, 
	const double *vl,
    double** DM_R)
{
	ModuleBase::TITLE("Gint_k","cal_force_k");
	ModuleBase::timer::tick("Gint_k","cal_force_k");

	int nnrg = GlobalC::GridT.nnrg;

	if(GlobalV::OUT_LEVEL != "m") GlobalV::ofs_running << " LNNR.nnrg in cal_force_k = " << nnrg << std::endl;
	assert(nnrg>=0);

	const double delta_r = GlobalC::ORB.dr_uniform;
	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	const int max_size = GlobalC::GridT.max_atom;
	// how many meshcells in bigcell.
	const int bxyz = GlobalC::GridT.bxyz;

	if(max_size!=0)
	{
		assert(this->ncxyz!=0);
		const double dv = GlobalC::ucell.omega/this->ncxyz;
		const int ncyz = GlobalC::pw.ncy*GlobalC::pw.nczp; 

		for(int i=0; i<nbx; i++)
		{
			const int ibx = i*GlobalC::pw.bx; // mohan add 2012-03-25
			for(int j=0; j<nby; j++)
			{
				const int jby = j*GlobalC::pw.by; // mohan add 2012-03-25
				for(int k=nbz_start; k<nbz_start+nbz; k++)
				{
					const int kbz = k*GlobalC::pw.bz-GlobalC::pw.nczp_start; //mohan add 2012-03-25
					const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
					const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
					if(na_grid==0) continue;

					int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);

                    int * block_iw, * block_index, * block_size;
                    Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);
					bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);
					//---------------------------------
					// get the wave functions in this
					// grid.
					//---------------------------------

                    // set up band matrix psir_ylm and psir_DM
                    const int LD_pool = max_size*GlobalC::ucell.nwmax;

                    Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_x(GlobalC::pw.bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_y(GlobalC::pw.bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_z(GlobalC::pw.bxyz, LD_pool);

                    Gint_Tools::cal_dpsir_ylm(
                        na_grid, grid_index, delta_r,
                        block_index, block_size, 
                        cal_flag,
                        psir_ylm.ptr_2D,
                        dpsir_ylm_x.ptr_2D,
                        dpsir_ylm_y.ptr_2D,
                        dpsir_ylm_z.ptr_2D
                    );

					double *vldr3 = Gint_Tools::get_vldr3(vl, ncyz, ibx, jby, kbz, dv);
					const Gint_Tools::Array_Pool<double> psir_vlbr3 
						= Gint_Tools::get_psir_vlbr3(na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);
					const Gint_Tools::Array_Pool<double> psir_vlbr3_DMR
						= Gint_Tools::get_psir_vlbr3_DMR(
							grid_index, na_grid,
							block_index, block_size,
							cal_flag, psir_vlbr3.ptr_2D,
							dpsir_ylm_x.ptr_2D,
							dpsir_ylm_y.ptr_2D,
							dpsir_ylm_z.ptr_2D,
							GlobalC::GridT, DM_R);

					if(isforce)
					{
                        this-> cal_meshball_force(
							grid_index, na_grid, 
                            block_size, block_index,
                            psir_vlbr3_DMR.ptr_2D, 
                            dpsir_ylm_x.ptr_2D, 
                            dpsir_ylm_y.ptr_2D, 
                            dpsir_ylm_z.ptr_2D, 
                            fvl_dphi);
					}
					if(isstress)
					{
						Gint_Tools::Array_Pool<double> dpsir_ylm_xx(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_xy(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_xz(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_yy(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_yz(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> dpsir_ylm_zz(GlobalC::pw.bxyz, LD_pool);
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
                            psir_vlbr3_DMR.ptr_2D, 
                            dpsir_ylm_xx.ptr_2D, 
                            dpsir_ylm_xy.ptr_2D, 
                            dpsir_ylm_xz.ptr_2D,
                            dpsir_ylm_yy.ptr_2D, 
                            dpsir_ylm_yz.ptr_2D, 
                            dpsir_ylm_zz.ptr_2D,
                            svl_dphi);
					}
					free(vindex);			vindex=nullptr;
					delete[] block_iw;
					delete[] block_index;
					delete[] block_size;
					for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
						free(cal_flag[ib]);
					free(cal_flag);			cal_flag=nullptr;
				}// int k
			}// int j
		} // int i
	}

	ModuleBase::timer::tick("Gint_k","cal_force_k");
	return;
}

void Gint_k::cal_meshball_force(
    const int grid_index,
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    ModuleBase::matrix &force)
{

	const int inc=1;
    for(int ia1=0;ia1<na_grid;ia1++)
    {
        const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + ia1;
        const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom

        for(int ib=0;ib<GlobalC::pw.bxyz;ib++)
        {
            const double rx = ddot_(&block_size[ia1], &psir_vlbr3_DMR[ib][block_index[ia1]], &inc, &dpsir_x[ib][block_index[ia1]], &inc);
            force(iat,0)+=rx*2.0;
            const double ry = ddot_(&block_size[ia1], &psir_vlbr3_DMR[ib][block_index[ia1]], &inc, &dpsir_y[ib][block_index[ia1]], &inc);
            force(iat,1)+=ry*2.0;
            const double rz = ddot_(&block_size[ia1], &psir_vlbr3_DMR[ib][block_index[ia1]], &inc, &dpsir_z[ib][block_index[ia1]], &inc);
            force(iat,2)+=rz*2.0;
          
        }
    }
	
	return;
}

void Gint_k::cal_meshball_stress(
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DMR,
    const double*const*const dpsir_xx,
    const double*const*const dpsir_xy,
    const double*const*const dpsir_xz,
    const double*const*const dpsir_yy,
    const double*const*const dpsir_yz,
    const double*const*const dpsir_zz,
    ModuleBase::matrix &stress)
{
    constexpr int inc=1;
    for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    {
        const double rxx = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_xx[ib], &inc);
        stress(0,0)+=rxx*2.0;
        const double rxy = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_xy[ib], &inc);
        stress(0,1)+=rxy*2.0;
        const double rxz = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_xz[ib], &inc);
        stress(0,2)+=rxz*2.0;
        const double ryy = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_yy[ib], &inc);
        stress(1,1)+=ryy*2.0;
        const double ryz = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_yz[ib], &inc);
        stress(1,2)+=ryz*2.0;
        const double rzz = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_zz[ib], &inc);
        stress(2,2)+=rzz*2.0;
    }
    return;
}