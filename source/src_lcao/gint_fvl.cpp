#include "gint_k.h"
#include "../src_pw/global.h"
#include "global_fp.h" // mohan add 2021-01-30

#include "../module_base/ylm.h"
#include "../module_base/timer.h"

void Gint::gint_kernel_force(
	const int na_grid,
	const int grid_index,
	const double delta_r,
	double* vldr3,
	const int LD_pool,
	double** DM_in,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix* fvl_dphi,
    ModuleBase::matrix* svl_dphi)
{
    //prepare block information
	int * block_iw, * block_index, * block_size;
	Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);
	bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

    //evaluate psi and dpsi on grids
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

    //calculating f_mu(r) = v(r)*psi_mu(r)*dv
	const Gint_Tools::Array_Pool<double> psir_vlbr3 
		= Gint_Tools::get_psir_vlbr3(na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);

	Gint_Tools::Array_Pool<double> psir_vlbr3_DM(GlobalC::pw.bxyz, LD_pool);
	ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_DM.ptr_1D, GlobalC::pw.bxyz*LD_pool);

	//calculating g_mu(r) = sum_nu rho_mu,nu f_nu(r)
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		Gint_Tools::mult_psi_DM(
			na_grid, LD_pool,
			block_iw, block_size,
			block_index, cal_flag,
			psir_vlbr3.ptr_2D,
			psir_vlbr3_DM.ptr_2D,
			DM_in, 2);
	}
	else
	{
		Gint_Tools::mult_psi_DMR(
			grid_index, na_grid,
			block_index, block_size,
			cal_flag, GlobalC::GridT,
			psir_vlbr3.ptr_2D,
			psir_vlbr3_DM.ptr_2D,
			DM_in[GlobalV::CURRENT_SPIN], 2);
	}

	if(isforce)
	{
        //do integration to get force
		this-> cal_meshball_force(
			grid_index, na_grid, 
			block_size, block_index,
			psir_vlbr3_DM.ptr_2D, 
			dpsir_ylm_x.ptr_2D, 
			dpsir_ylm_y.ptr_2D, 
			dpsir_ylm_z.ptr_2D, 
			fvl_dphi);
	}
	if(isstress)
	{
        //calculating g_mu(r)*(r-R) where R is the location of atom
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
        //do integration to get stress
		this-> cal_meshball_stress(na_grid, block_index,
			psir_vlbr3_DM.ptr_2D, 
			dpsir_ylm_xx.ptr_2D, 
			dpsir_ylm_xy.ptr_2D, 
			dpsir_ylm_xz.ptr_2D,
			dpsir_ylm_yy.ptr_2D, 
			dpsir_ylm_yz.ptr_2D, 
			dpsir_ylm_zz.ptr_2D,
			svl_dphi);
	}

    //release memories
	delete[] block_iw;
	delete[] block_index;
	delete[] block_size;
	for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
	{
		delete[] cal_flag[ib];
	}
	delete[] cal_flag;
}

void Gint::cal_meshball_force(
    const int grid_index,
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    ModuleBase::matrix *force)
{

	const int inc=1;
    for(int ia1=0;ia1<na_grid;ia1++)
    {
        const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + ia1;
        const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom

        for(int ib=0;ib<GlobalC::pw.bxyz;ib++)
        {
            const double rx = ddot_(&block_size[ia1], &psir_vlbr3_DMR[ib][block_index[ia1]], &inc, &dpsir_x[ib][block_index[ia1]], &inc);
            force[0](iat,0)+=rx*2.0;
            const double ry = ddot_(&block_size[ia1], &psir_vlbr3_DMR[ib][block_index[ia1]], &inc, &dpsir_y[ib][block_index[ia1]], &inc);
            force[0](iat,1)+=ry*2.0;
            const double rz = ddot_(&block_size[ia1], &psir_vlbr3_DMR[ib][block_index[ia1]], &inc, &dpsir_z[ib][block_index[ia1]], &inc);
            force[0](iat,2)+=rz*2.0;
          
        }
    }
	
	return;
}

void Gint::cal_meshball_stress(
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DMR,
    const double*const*const dpsir_xx,
    const double*const*const dpsir_xy,
    const double*const*const dpsir_xz,
    const double*const*const dpsir_yy,
    const double*const*const dpsir_yz,
    const double*const*const dpsir_zz,
    ModuleBase::matrix *stress)
{
    constexpr int inc=1;
    for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    {
        const double rxx = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_xx[ib], &inc);
        stress[0](0,0)+=rxx*2.0;
        const double rxy = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_xy[ib], &inc);
        stress[0](0,1)+=rxy*2.0;
        const double rxz = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_xz[ib], &inc);
        stress[0](0,2)+=rxz*2.0;
        const double ryy = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_yy[ib], &inc);
        stress[0](1,1)+=ryy*2.0;
        const double ryz = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_yz[ib], &inc);
        stress[0](1,2)+=ryz*2.0;
        const double rzz = ddot_(&block_index[na_grid], psir_vlbr3_DMR[ib], &inc, dpsir_zz[ib], &inc);
        stress[0](2,2)+=rzz*2.0;
    }
    return;
}