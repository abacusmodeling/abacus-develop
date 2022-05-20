#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "gint_tools.h"

void Gint_k::gint_kernel_rho(
	const int na_grid,
	const int grid_index,
	const double delta_r,
	int* vindex,
	const int LD_pool,
	Gint_inout *inout)
{
	//prepare block information
	int * block_iw, * block_index, * block_size;
	Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);
	bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

	//evaluate psi on grids
	Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);
	Gint_Tools::cal_psir_ylm(
		na_grid, grid_index, delta_r,
		block_index, block_size, 
		cal_flag,
		psir_ylm.ptr_2D);

	//setting variables
	double** DM_R = inout->DM_R;
	Charge* chr = inout->chr;

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		//calculating g_mu(r) = sum_nu rho_mu,nu psi_nu(r)
		const Gint_Tools::Array_Pool<double> psir_DMR
			= Gint_Tools::mult_psi_DMR(
				grid_index, na_grid,
				block_index, block_size,
				cal_flag, GlobalC::GridT,
				psir_ylm.ptr_2D,
				DM_R[is], 1);

		//do sum_mu g_mu(r)psi_mu(r) to get electron density on grid
		this->cal_meshball_rho(
			na_grid,
			block_index,
			vindex, 
			psir_ylm.ptr_2D,
			psir_DMR.ptr_2D,
			chr->rho[is]);
	}
	delete[] block_iw;
	delete[] block_index;
	delete[] block_size;
	for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
	{
		delete[] cal_flag[ib];
	}
	delete[] cal_flag;
}

void Gint_k::cal_meshball_rho(
	const int na_grid,
	int* block_index,
	int* vindex,
	double** psir_ylm,
	double** psir_DMR,
	double* rho)
{		
	const int inc = 1;
	// sum over mu to get density on grid
	for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
	{
		double r=ddot_(&block_index[na_grid], psir_ylm[ib], &inc, psir_DMR[ib], &inc);
		const int grid = vindex[ib];
		rho[ grid ] += r;
	}
}


