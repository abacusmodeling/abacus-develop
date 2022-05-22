//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

// can be done by GPU
void Gint_Gamma::gint_kernel_rho(
	const int na_grid,
	const int grid_index,
	const double delta_r,
	int* vindex,
	const int LD_pool,
	Gint_inout *inout) const		// rho[GlobalV::NSPIN][GlobalC::pw.nrxx]
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

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        Gint_Tools::Array_Pool<double> psir_DM
			= Gint_Tools::mult_psi_DM(
				na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				psir_ylm.ptr_2D,
				inout->DM[is], 1);

		this->cal_meshball_rho(na_grid, block_index,
			vindex, psir_ylm.ptr_2D,
			psir_DM.ptr_2D, inout->chr->rho[is]);
    } // end is
	delete[] block_iw;
	delete[] block_index;
	delete[] block_size;

	for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
	{
		delete[] cal_flag[ib];
	}
	delete[] cal_flag;

}

// calculate charge density
void Gint_Gamma::cal_rho(Gint_inout *inout)
{
    ModuleBase::TITLE("Gint_Gamma","cal_rho");
    ModuleBase::timer::tick("Gint_Gamma","cal_rho");

    max_size = GlobalC::GridT.max_atom;
	const int LD_pool = max_size*GlobalC::ucell.nwmax;

	if(max_size)
    {
#ifdef __MKL
   		const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(std::max(1,mkl_threads/GlobalC::GridT.nbx));		// Peize Lin update 2021.01.20
#endif
		
#ifdef _OPENMP
		#pragma omp parallel
#endif
		{		
			const int nbx = GlobalC::GridT.nbx;
			const int nby = GlobalC::GridT.nby;
			const int nbz_start = GlobalC::GridT.nbzp_start;
			const int nbz = GlobalC::GridT.nbzp;
		
			const int ncyz = GlobalC::pw.ncy*GlobalC::pw.nczp; // mohan add 2012-03-25
            
            // it's a uniform grid to save orbital values, so the delta_r is a constant.
            const double delta_r = GlobalC::ORB.dr_uniform;		
#ifdef _OPENMP
			#pragma omp for
#endif
			for (int i=0; i<nbx; i++)
			{
				const int ibx = i*GlobalC::pw.bx;
				for (int j=0; j<nby; j++)
				{
					const int jby = j*GlobalC::pw.by;
					for (int k=nbz_start; k<nbz_start+nbz; k++)
					{
						const int kbz = k*GlobalC::pw.bz-GlobalC::pw.nczp_start;
		
						const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
		
						// get the value: how many atoms has orbital value on this grid.
						const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
						if(na_grid==0) continue;				
						
						// here vindex refers to local potentials
						if(inout->job == Gint_Tools::job_type::rho)
						{
							int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);						
							this->gint_kernel_rho(na_grid, grid_index,
								delta_r, vindex, LD_pool, 
								inout); 
							delete[] vindex;
						}
					}// k
				}// j
			}// i
		} // end of #pragma omp parallel
			
#ifdef __MKL
   		mkl_set_num_threads(mkl_threads);
#endif
    } // end of if(max_size)

    ModuleBase::timer::tick("Gint_Gamma","cal_rho");
    return;
}

void Gint_Gamma::cal_meshball_rho(
	const int na_grid,
	const int*const block_index,
	const int*const vindex,
	const double*const*const psir_ylm,
	double** psir_DM,
	double* rho) const
{		
	const int inc = 1;
	// sum over mu to get density on grid
	for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
	{
		double r=ddot_(&block_index[na_grid], psir_ylm[ib], &inc, psir_DM[ib], &inc);
		const int grid = vindex[ib];
		rho[ grid ] += r;
	}
}