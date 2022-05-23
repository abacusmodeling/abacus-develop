#include "gint_gamma.h"
#include "../src_pw/global.h"
#include "../module_base/ylm.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_base/timer.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint_Gamma::Gint_Gamma()
{
   
    sender_index_size = 1;
	sender_local_index = new int[1];
    sender_size_process = new int[1];
    sender_displacement_process = new int[1];
    sender_size=1;
    sender_buffer=new double[1];

    receiver_index_size=1;
    receiver_global_index = new int[1];
    receiver_size_process = new int[1];
    receiver_displacement_process = new int[1];
    receiver_size=1;
    receiver_buffer=new double[1];
}

Gint_Gamma::~Gint_Gamma()
{

    delete[] sender_local_index;
    delete[] sender_size_process;
    delete[] sender_displacement_process;
    delete[] sender_buffer;

    delete[] receiver_global_index;
    delete[] receiver_size_process;
    delete[] receiver_displacement_process;
    delete[] receiver_buffer;
}

// calculate charge density
void Gint_Gamma::cal_gint_gamma(Gint_inout *inout)
{
    ModuleBase::TITLE("Gint_Gamma","cal_gint_gamma");
    ModuleBase::timer::tick("Gint_Gamma","cal_gint_gamma");

    const int max_size = GlobalC::GridT.max_atom;
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
			const double dv = GlobalC::ucell.omega/GlobalC::pw.ncxyz;
            
            // it's a uniform grid to save orbital values, so the delta_r is a constant.
            const double delta_r = GlobalC::ORB.dr_uniform;		
#ifdef _OPENMP
			ModuleBase::matrix fvl_dphi_thread;
			ModuleBase::matrix svl_dphi_thread;
			if(inout->job==Gint_Tools::job_type::force)
			{
				if(inout->isforce)
				{
					fvl_dphi_thread.create(inout->fvl_dphi->nr,inout->fvl_dphi->nc);
					fvl_dphi_thread.zero_out();
				}
				if(inout->isstress)
				{
					svl_dphi_thread.create(inout->svl_dphi->nr,inout->svl_dphi->nc);
					svl_dphi_thread.zero_out();
				}
			}

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
						if(inout->job == Gint_Tools::job_type::force)
						{
							double* vldr3 = Gint_Tools::get_vldr3(inout->vl, ncyz, ibx, jby, kbz, dv);
                            #ifdef _OPENMP
                                this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool,
                                    inout->DM[GlobalV::CURRENT_SPIN],inout->isforce,
                                    inout->isstress, &fvl_dphi_thread, &svl_dphi_thread);
                            #else
                                this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool,
                                    inout->DM[GlobalV::CURRENT_SPIN],inout->isforce,
                                    inout->isstress, inout->fvl_dphi, inout->svl_dphi);
                            #endif
							delete [] vldr3;
						}
					}// k
				}// j
			}// i

#ifdef _OPENMP
			#pragma omp critical(gint_k)
			if(inout->job==Gint_Tools::job_type::force)
			{
				if(inout->isforce)
				{
					inout->fvl_dphi[0]+=fvl_dphi_thread;
				}
				if(inout->isstress)
				{
					inout->svl_dphi[0]+=svl_dphi_thread;
				}
			}
#endif

		} // end of #pragma omp parallel
			
#ifdef __MKL
   		mkl_set_num_threads(mkl_threads);
#endif
    } // end of if(max_size)

    ModuleBase::timer::tick("Gint_Gamma","cal_gint_gamma");
    return;
}