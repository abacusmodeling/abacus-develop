#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "../module_base/memory.h"
#include "../module_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint_k::Gint_k()
{
    pvpR_alloc_flag = false;
    spin_now = -1; // for a start value, must not equal 1,2 or 4.
}

Gint_k::~Gint_k()
{
    
}

void Gint_k::prep_grid(
	const int &nbx_in,
	const int &nby_in,
	const int &nbz_in,
	const int &nbz_start_in,
    const int& ncxyz_in)
{
	ModuleBase::TITLE(GlobalV::ofs_running,"Gint_k","prep_grid");
    
    this->nbx = nbx_in;
	this->nby = nby_in;
	this->nbz = nbz_in;
	this->ncxyz = ncxyz_in;
	this->nbz_start = nbz_start_in;
	assert(nbx>0);
	assert(nby>0);
	assert(nbz>0);
	assert(ncxyz>0);

	assert( GlobalC::ucell.omega > 0.0);

	return;
}

void Gint_k::cal_gint_k(Gint_inout *inout)
{
	ModuleBase::TITLE("Gint_k","cal_gint_k");
    ModuleBase::timer::tick("Gint_k", "cal_gint_k");

	const int max_size = GlobalC::GridT.max_atom;

    if(max_size!=0)
    {
#ifdef __MKL
		const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(1);
#endif

#ifdef _OPENMP
    	#pragma omp parallel
#endif
		{
			const int ncyz = GlobalC::pw.ncy*GlobalC::pw.nczp; // mohan add 2012-03-25
			const double dv = GlobalC::ucell.omega/this->ncxyz;
			
			// it's a uniform grid to save orbital values, so the delta_r is a constant.
			const double delta_r = GlobalC::ORB.dr_uniform;	
#ifdef _OPENMP
			double* pvpR_reduced_thread;
			if(inout->job==Gint_Tools::job_type::vlocal)
			{
        		pvpR_reduced_thread = new double[GlobalC::GridT.nnrg];
        		ModuleBase::GlobalFunc::ZEROS(pvpR_reduced_thread, GlobalC::GridT.nnrg);
			}

    		#pragma omp for
#endif
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

						// get the value: how many atoms has orbital value on this grid.
						const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
						if(na_grid==0) continue;				

						// set up band matrix psir_ylm and psir_DM
						const int LD_pool = max_size*GlobalC::ucell.nwmax;
						
						switch(inout->job)
						{
							case Gint_Tools::job_type::rho:
								int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);
								this->gint_kernel_rho(na_grid, grid_index, delta_r, vindex, LD_pool, inout);
								delete[] vindex;
								break;
							case Gint_Tools::job_type::force:
								double *vldr3 = Gint_Tools::get_vldr3(inout->vl, ncyz, ibx, jby, kbz, dv);
								this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool, inout);
								delete[] vldr3;
								break;
						}

					}// int k
				}// int j
			} // int i
		} // end of #pragma omp parallel
#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
    } // end of if (max_size)	

	ModuleBase::timer::tick("Gint_k","cal_gint_k");
	return;
}
