//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_tools.h"
#include "../src_pw/global.h"
#include "global_fp.h"
#include "../module_base/ylm.h"
#include <cmath>

namespace Gint_Tools
{
	// here vindex refers to local potentials
	int* get_vindex(
		const int ncyz,
		const int ibx,
		const int jby,
		const int kbz)
	{
		int *vindex = (int*)malloc(GlobalC::pw.bxyz*sizeof(int));
		int bindex=0;
		// z is the fastest,
		// ipart can be obtained by using a previously stored array
		for(int ii=0; ii<GlobalC::pw.bx; ii++)
		{
			const int ipart = (ibx+ii)*ncyz;
			for(int jj=0; jj<GlobalC::pw.by; jj++)
			{
				// jpart can be obtained by using a previously stored array
				const int jpart = (jby+jj)*GlobalC::pw.nczp + ipart;
				for(int kk=0; kk<GlobalC::pw.bz; kk++)
				{
					vindex[bindex] = kbz+kk + jpart;
					++bindex;
				}
			}
		}
		return vindex;
	}

	// extract the local potentials.
	double* get_vldr3(
		const double*const vlocal,		// vlocal[ir]
		const int ncyz,
		const int ibx,
		const int jby,
		const int kbz,
		const double dv)
	{
		// set the index for obtaining local potentials
		int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);	
		double *vldr3 = (double*)malloc(GlobalC::pw.bxyz*sizeof(double));					
		for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
		{
			vldr3[ib]=vlocal[vindex[ib]] * dv;
		}
		free(vindex);	vindex=nullptr;
		return vldr3;
	}

	void get_block_info(
		const int na_grid,
		const int grid_index,
		int * &block_iw,
		int * &block_index,
		int * &block_size	
	)
	{
		block_iw = new int[na_grid];
		block_index = new int[na_grid+1];
		block_size = new int[na_grid];

		block_index[0] = 0;
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + id;
			const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom
			const int it=GlobalC::ucell.iat2it[ iat ]; // index of atom type
			const int ia=GlobalC::ucell.iat2ia[ iat ]; // index of atoms within each type
			const int start=GlobalC::ucell.itiaiw2iwt(it, ia, 0); // the index of the first wave function for atom (it,ia)
			block_iw[id]=GlobalC::GridT.trace_lo[start];
			block_index[id+1] = block_index[id]+GlobalC::ucell.atoms[it].nw;
			block_size[id]=GlobalC::ucell.atoms[it].nw;	
		}
	}

	void get_block_info(
		const int na_grid,
		const int grid_index,
		int * &block_iw,
		int * &block_index,
		int * &block_size,
		int * &at
	)
	{
		block_iw = new int[na_grid];
		block_index = new int[na_grid+1];
		block_size = new int[na_grid];
		at = new int[na_grid];

		block_index[0] = 0;
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + id;
			const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom
			const int it=GlobalC::ucell.iat2it[ iat ]; // index of atom type
			const int ia=GlobalC::ucell.iat2ia[ iat ]; // index of atoms within each type
			const int start=GlobalC::ucell.itiaiw2iwt(it, ia, 0); // the index of the first wave function for atom (it,ia)
			block_iw[id]=GlobalC::GridT.trace_lo[start];
			block_index[id+1] = block_index[id]+GlobalC::ucell.atoms[it].nw;
			block_size[id]=GlobalC::ucell.atoms[it].nw;	
			at[id] = iat;
		}
	}
	// whether the atom-grid distance is larger than cutoff
	bool** get_cal_flag(
		const int na_grid, 			// number of atoms on this grid 
		const int grid_index)		// 1d index of FFT index (i,j,k) 
	{
		bool** cal_flag = (bool**)malloc(GlobalC::pw.bxyz*sizeof(bool*));
		for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
			cal_flag[ib] = (bool*)malloc(na_grid*sizeof(bool));

		for (int id=0; id<na_grid; id++)
		{
			// there are two parameters we want to know here:
			// in which bigcell of the meshball the atom in?
			// what's the cartesian coordinate of the bigcell?
			const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + id;
			const int iat=GlobalC::GridT.which_atom[mcell_index];		
			const int it=GlobalC::ucell.iat2it[iat];

			// meshball_positions should be the bigcell position in meshball
			// to the center of meshball.
			// calculated in cartesian coordinates
			// the std::vector from the grid which is now being operated to the atom position.
			// in meshball language, is the std::vector from imcell to the center cel, plus
			// tau_in_bigcell.
			const int imcell=GlobalC::GridT.which_bigcell[mcell_index];
			const double mt[3] = {
				GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0],
				GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1],
				GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
			{
				// meshcell_pos: z is the fastest
				const double dr[3] = {
					GlobalC::GridT.meshcell_pos[ib][0] + mt[0],
					GlobalC::GridT.meshcell_pos[ib][1] + mt[1],
					GlobalC::GridT.meshcell_pos[ib][2] + mt[2]};
				const double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);	// distance between atom and grid

				if(distance > (GlobalC::ORB.Phi[it].getRcut()-1.0e-15))
					cal_flag[ib][id]=false;
				else
					cal_flag[ib][id]=true;
			}// end ib
		}// end id
		return cal_flag;
	}

	void cal_psir_ylm(
		const int na_grid, 					// number of atoms on this grid 
		const int grid_index, 				// 1d index of FFT index (i,j,k) 
		const double delta_r, 				// delta_r of the uniform FFT grid
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,
		double*const*const psir_ylm) 	// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	{
		for (int id=0; id<na_grid; id++)
		{
			// there are two parameters we want to know here:
			// in which bigcell of the meshball the atom is in?
			// what's the cartesian coordinate of the bigcell?
			const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + id;

			const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom
			const int it=GlobalC::ucell.iat2it[iat]; // index of atom type
			const Atom*const atom=&GlobalC::ucell.atoms[it];

			// meshball_positions should be the bigcell position in meshball
			// to the center of meshball.
			// calculated in cartesian coordinates
			// the std::vector from the grid which is now being operated to the atom position.
			// in meshball language, is the std::vector from imcell to the center cel, plus
			// tau_in_bigcell.
			const int imcell=GlobalC::GridT.which_bigcell[mcell_index];
			const double mt[3] = {
				GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0],
				GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1],
				GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2]};

			// number of grids in each big cell (bxyz)
			for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
			{
				double *p=&psir_ylm[ib][block_index[id]];
				if(!cal_flag[ib][id]) 
				{
					ModuleBase::GlobalFunc::ZEROS(p, block_size[id]);
				}
				else
				{
					// meshcell_pos: z is the fastest
					const double dr[3] = {
						GlobalC::GridT.meshcell_pos[ib][0] + mt[0],
						GlobalC::GridT.meshcell_pos[ib][1] + mt[1],
						GlobalC::GridT.meshcell_pos[ib][2] + mt[2]};	
					double distance = std::sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );	// distance between atom and grid
					//if(distance[id] > GlobalC::GridT.orbital_rmax) continue;
					if (distance < 1.0E-9) distance += 1.0E-9;
					
					//------------------------------------------------------
					// spherical harmonic functions Ylm
					//------------------------------------------------------
					std::vector<double> ylma;
					//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
					ModuleBase::Ylm::sph_harm ( GlobalC::ucell.atoms[it].nwl,
							dr[0] / distance,
							dr[1] / distance,
							dr[2] / distance,
							ylma);
					// these parameters are related to interpolation
					// because once the distance from atom to grid point is known,
					// we can obtain the parameters for interpolation and
					// store them first! these operations can save lots of efforts.
					const double position = distance / delta_r;
					const int ip = static_cast<int>(position);
					const double dx = position - ip;
					const double dx2 = dx * dx;
					const double dx3 = dx2 * dx;

					const double c3 = 3.0*dx2-2.0*dx3;
					const double c1 = 1.0-c3;
					const double c2 = (dx-2.0*dx2+dx3)*delta_r;
					const double c4 = (dx3-dx2)*delta_r;

					double phi=0;
					for (int iw=0; iw< atom->nw; ++iw, ++p)
					{
						if ( atom->iw2_new[iw] )
						{
							const Numerical_Orbital_Lm &philn = GlobalC::ORB.Phi[it].PhiLN(
									atom->iw2l[iw],
									atom->iw2n[iw]);
							phi = c1*philn.psi_uniform[ip] + c2*philn.dpsi_uniform[ip]			 // radial wave functions
								+ c3*philn.psi_uniform[ip+1] + c4*philn.dpsi_uniform[ip+1];
						}
						*p=phi * ylma[atom->iw2_ylm[iw]];
					} // end iw
				}// end distance<=(GlobalC::ORB.Phi[it].getRcut()-1.0e-15)
			}// end ib
		}// end id
		return;
	}

	void cal_dpsir_ylm(
		const int na_grid, 					// number of atoms on this grid 
		const int grid_index, 				// 1d index of FFT index (i,j,k) 
		const double delta_r, 				// delta_r of the uniform FFT grid
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,    // cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		double*const*const psir_ylm,
		double*const*const dpsir_ylm_x,
		double*const*const dpsir_ylm_y,
		double*const*const dpsir_ylm_z)
	{
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index = GlobalC::GridT.bcell_start[grid_index] + id;
			const int imcell = GlobalC::GridT.which_bigcell[mcell_index];
			int iat = GlobalC::GridT.which_atom[mcell_index];
			const int it = GlobalC::ucell.iat2it[iat];
			const int ia = GlobalC::ucell.iat2ia[iat];
			Atom *atom = &GlobalC::ucell.atoms[it];

			const double mt[3]={
				GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0],
				GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1],
				GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
			{
				double*const p_psi=&psir_ylm[ib][block_index[id]];
				double*const p_dpsi_x=&dpsir_ylm_x[ib][block_index[id]];
				double*const p_dpsi_y=&dpsir_ylm_y[ib][block_index[id]];
				double*const p_dpsi_z=&dpsir_ylm_z[ib][block_index[id]];
				if(!cal_flag[ib][id]) 
				{
					ModuleBase::GlobalFunc::ZEROS(p_psi, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_x, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_y, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_z, block_size[id]);
				}
				else
				{
					const double dr[3]={						// vectors between atom and grid
						GlobalC::GridT.meshcell_pos[ib][0] + mt[0],
						GlobalC::GridT.meshcell_pos[ib][1] + mt[1],
						GlobalC::GridT.meshcell_pos[ib][2] + mt[2]};
					double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

					//array to store spherical harmonics and its derivatives
					std::vector<double> rly;
					std::vector<std::vector<double>> grly;
					ModuleBase::Ylm::grad_rl_sph_harm(GlobalC::ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly);
					if(distance < 1e-9)  distance = 1e-9;

					const double position = distance / delta_r;
							
					const double iq = static_cast<int>(position);
					const double x0 = position - iq;
					const double x1 = 1.0 - x0;
					const double x2 = 2.0 - x0;
					const double x3 = 3.0 - x0;
					const double x12 = x1*x2 / 6;
					const double x03 = x0*x3 / 2;
					
					double tmp, dtmp;

					for (int iw=0; iw< atom->nw; ++iw)
					{
						// this is a new 'l', we need 1D orbital wave
						// function from interpolation method.
						if ( atom->iw2_new[iw] )
						{
							const Numerical_Orbital_Lm &philn = GlobalC::ORB.Phi[it].PhiLN(
									atom->iw2l[iw],
									atom->iw2n[iw]);

							//if ( iq[id] >= philn.nr_uniform-4)
							if ( iq >= philn.nr_uniform-4)
							{
								tmp = dtmp = 0.0;
							}
							else
							{
								// use Polynomia Interpolation method to get the 
								// wave functions

								tmp = x12*(philn.psi_uniform[iq]*x3
										+philn.psi_uniform[iq+3]*x0)
									+ x03*(philn.psi_uniform[iq+1]*x2
											-philn.psi_uniform[iq+2]*x1);

								dtmp = x12*(philn.dpsi_uniform[iq]*x3
										+philn.dpsi_uniform[iq+3]*x0)
										+ x03*(philn.dpsi_uniform[iq+1]*x2
											-philn.dpsi_uniform[iq+2]*x1);
							}
						}//new l is used.
						
						// get the 'l' of this localized wave function
						const int ll = atom->iw2l[iw];
						const int idx_lm = atom->iw2_ylm[iw];

						const double rl = pow(distance, ll);

						// 3D wave functions
						p_psi[iw] = tmp * rly[idx_lm] / rl;

						// derivative of wave functions with respect to atom positions.
						const double tmpdphi_rly = (dtmp  - tmp * ll / distance) / rl * rly[idx_lm] / distance;
						const double tmprl = tmp/rl;

						p_dpsi_x[iw] = tmpdphi_rly * dr[0]  + tmprl * grly[idx_lm][0];
						p_dpsi_y[iw] = tmpdphi_rly * dr[1]  + tmprl * grly[idx_lm][1];
						p_dpsi_z[iw] = tmpdphi_rly * dr[2]  + tmprl * grly[idx_lm][2];
					}//iw
				}//else
			}	
		}

		return;
	}

	void cal_dpsirr_ylm(
		const int na_grid, 					// number of atoms on this grid 
		const int grid_index, 				// 1d index of FFT index (i,j,k) 
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,    // cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		double*const*const dpsir_ylm_x,
		double*const*const dpsir_ylm_y,
		double*const*const dpsir_ylm_z,
		double*const*const dpsir_ylm_xx,
		double*const*const dpsir_ylm_xy,
		double*const*const dpsir_ylm_xz,
		double*const*const dpsir_ylm_yy,
		double*const*const dpsir_ylm_yz,
		double*const*const dpsir_ylm_zz)
	{
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index = GlobalC::GridT.bcell_start[grid_index] + id;
			const int imcell = GlobalC::GridT.which_bigcell[mcell_index];
			int iat = GlobalC::GridT.which_atom[mcell_index];
			const int it = GlobalC::ucell.iat2it[iat];
			Atom *atom = &GlobalC::ucell.atoms[it];

			const double mt[3]={
				GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0],
				GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1],
				GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
			{
				double*const p_dpsi_x=&dpsir_ylm_x[ib][block_index[id]];
				double*const p_dpsi_y=&dpsir_ylm_y[ib][block_index[id]];
				double*const p_dpsi_z=&dpsir_ylm_z[ib][block_index[id]];
				double*const p_dpsi_xx=&dpsir_ylm_xx[ib][block_index[id]];
				double*const p_dpsi_xy=&dpsir_ylm_xy[ib][block_index[id]];
				double*const p_dpsi_xz=&dpsir_ylm_xz[ib][block_index[id]];
				double*const p_dpsi_yy=&dpsir_ylm_yy[ib][block_index[id]];
				double*const p_dpsi_yz=&dpsir_ylm_yz[ib][block_index[id]];
				double*const p_dpsi_zz=&dpsir_ylm_zz[ib][block_index[id]];
				if(!cal_flag[ib][id]) 
				{
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xx, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_yy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_yz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_zz, block_size[id]);
				}
				else
				{
					const double dr[3]={						// vectors between atom and grid
						GlobalC::GridT.meshcell_pos[ib][0] + mt[0],
						GlobalC::GridT.meshcell_pos[ib][1] + mt[1],
						GlobalC::GridT.meshcell_pos[ib][2] + mt[2]};

					for (int iw=0; iw< atom->nw; ++iw)
					{

						p_dpsi_xx[iw] = p_dpsi_x[iw]*dr[0];
						p_dpsi_xy[iw] = p_dpsi_x[iw]*dr[1];
						p_dpsi_xz[iw] = p_dpsi_x[iw]*dr[2];
						p_dpsi_yy[iw] = p_dpsi_y[iw]*dr[1];
						p_dpsi_yz[iw] = p_dpsi_y[iw]*dr[2];
						p_dpsi_zz[iw] = p_dpsi_z[iw]*dr[2];

					}//iw
				}//else
			}	
		}

		return;
	}

	// atomic basis sets
	// psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
	Gint_Tools::Array_Pool<double> get_psir_vlbr3(
		const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,	    	// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const vldr3,			    	// vldr3[GlobalC::pw.bxyz]
		const double*const*const psir_ylm)		    // psir_ylm[GlobalC::pw.bxyz][LD_pool]
	{
		Gint_Tools::Array_Pool<double> psir_vlbr3(GlobalC::pw.bxyz, LD_pool);
		for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
		{
			for(int ia=0; ia<na_grid; ++ia)
			{
				if(cal_flag[ib][ia])
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3.ptr_2D[ib][i]=psir_ylm[ib][i]*vldr3[ib];
					}
				}
				else
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3.ptr_2D[ib][i]=0;
					}
				}

			}
		}
		return psir_vlbr3;
	}

	Gint_Tools::Array_Pool<double> get_psir_vlbr3_DM(
		const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_iw,				    // block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,	    	// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const*const psir_vlbr3,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		const double*const*const DM)
	{
		constexpr char side='L', uplo='U';
		constexpr char transa='N', transb='N';
		constexpr double alpha_symm=1, alpha_gemm=1, beta=1;    
		constexpr int inc=1;

		Gint_Tools::Array_Pool<double> psir_vlbr3_DM(GlobalC::pw.bxyz, LD_pool);
		ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_DM.ptr_1D, GlobalC::pw.bxyz*LD_pool);

		for (int ia1=0; ia1<na_grid; ia1++)
		{
			const int iw1_lo=block_iw[ia1];
			for (int ia2=0; ia2<na_grid; ia2++)
			{
				int first_ib=0, last_ib=0;
				for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						first_ib=ib;
						break;
					}
				}
				for(int ib=GlobalC::pw.bxyz-1; ib>=0; --ib)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						last_ib=ib+1;
						break;
					}
				}
				const int ib_length=last_ib-first_ib;
				if(ib_length<=0) continue;

				int cal_pair_num=0;
				for(int ib=first_ib; ib<last_ib; ++ib)
				{
					cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
				}
				const int iw2_lo=block_iw[ia2];
				if(cal_pair_num>ib_length/4)
				{
					dgemm_(&transa, &transb, &block_size[ia2], &ib_length, &block_size[ia1], 
						&alpha_gemm, &DM[iw1_lo][iw2_lo], &GlobalC::GridT.lgd, 
						&psir_vlbr3[first_ib][block_index[ia1]], &LD_pool, 
						&beta, &psir_vlbr3_DM.ptr_2D[first_ib][block_index[ia2]], &LD_pool);
				}
				else
				{
					for(int ib=first_ib; ib<last_ib; ++ib)
					{
						if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
						{
							dgemv_(&transa, &block_size[ia2], &block_size[ia1], 
								&alpha_gemm, &DM[iw1_lo][iw2_lo], &GlobalC::GridT.lgd,
								&psir_vlbr3[ib][block_index[ia1]], &inc,
								&beta, &psir_vlbr3_DM.ptr_2D[ib][block_index[ia2]], &inc);
						}
					}
				}
			}// ia2       
		} // ia1  
		
		return psir_vlbr3_DM;
	}

//calculating (psi_DMR)_mu = sum_nu DMR_mu,nu psi_nu
//note : there is a difference between rho and force
//in calculating rho, due to symmetry, the summation over mu,nu
//can be done as sum_mu,mu + 2 sum_mu<nu, saving some time
//but for force, we cannot exchange the index
	Gint_Tools::Array_Pool<double> mult_psi_DMR(
		const int &grid_index,
		const int &na_grid,
		const int*const block_index,
		const int*const block_size,
		bool** cal_flag,
		const Grid_Technique &gt,
		double** psi, 
		double* DMR,
		const int job)
	{                       
		double *psi2, *psi2_dmr;
		int iwi, iww;
		const int LD_pool = GlobalC::GridT.max_atom*GlobalC::ucell.nwmax;
		Gint_Tools::Array_Pool<double> psi_DMR(GlobalC::pw.bxyz, LD_pool);
		ModuleBase::GlobalFunc::ZEROS(psi_DMR.ptr_1D, GlobalC::pw.bxyz*LD_pool);

		bool *all_out_of_range = new bool[na_grid];
		for(int ia=0; ia<na_grid; ++ia) //number of atoms
		{
			all_out_of_range[ia] = true;
			for(int ib=0; ib<gt.bxyz; ++ib) //number of small box in big box
			{
				if(cal_flag[ib][ia])
				{
					all_out_of_range[ia] = false;
					//break; //mohan add 2012-07-10
				}
			}
		}

		//parameters for lapack subroutiens
		const char trans='N';
		const double alpha=1.0, beta=1.0;
		const int inc=1;
		double alpha1;
		switch(job)
		{
			case 1:
				alpha1=2.0;
				break;
			case 2:
				alpha1=1.0;
				break;
			default:
				ModuleBase::WARNING_QUIT("psir_dmr","job can only be 1 or 2");
		}

		for (int ia1=0; ia1<na_grid; ia1++)
		{
			if(all_out_of_range[ia1]) continue;

			const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
			const int iat = gt.which_atom[mcell_index1];
			const int T1 = GlobalC::ucell.iat2it[iat];
			const int I1 = GlobalC::ucell.iat2ia[iat];
			Atom *atom1 = &GlobalC::ucell.atoms[T1];

			//~~~~~~~~~~~~~~~~
			// get cell R1.
			//~~~~~~~~~~~~~~~~
			const int id1 = gt.which_unitcell[mcell_index1];
			const int R1x = gt.ucell_index2x[id1];
			const int R1y = gt.ucell_index2y[id1];
			const int R1z = gt.ucell_index2z[id1];

			if(job==1) //density
			{
				const int idx1=block_index[ia1];
				int* find_start = GlobalC::GridT.find_R2[iat];
				int* find_end = GlobalC::GridT.find_R2[iat] + GlobalC::GridT.nad[iat];
				//ia2==ia1
				int cal_num=0;
				for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
				{
					if(cal_flag[ib][ia1])
					{
						++cal_num;
					}
				}

				int offset;
				if(cal_num>0)
				{
					//find offset				
					const int index = GlobalC::GridT.cal_RindexAtom(0, 0, 0, iat);
					offset = -1;
					for(int* find=find_start; find < find_end; find++)
					{
						//--------------------------------------------------------------
						// start positions of adjacent atom of 'iat'
						//--------------------------------------------------------------
						if( find[0] == index ) 
						{
							offset = find - find_start; // start positions of adjacent atom of 'iat'
							break;
						}
					}

					assert(offset!=-1);
					assert(offset < GlobalC::GridT.nad[iat]);				
				}

				if(cal_num>GlobalC::pw.bxyz/4)
				{				
					const int DM_start = GlobalC::GridT.nlocstartg[iat]+ GlobalC::GridT.find_R2st[iat][offset];					
					dgemm_(&trans, &trans, &block_size[ia1], &GlobalC::pw.bxyz, &block_size[ia1], &alpha,
						&DMR[DM_start], &block_size[ia1], 
						&psi[0][idx1], &LD_pool,  
						&beta, &psi_DMR.ptr_2D[0][idx1], &LD_pool);
				}
				else if(cal_num>0)
				{	
					const int DM_start = GlobalC::GridT.nlocstartg[iat]+ GlobalC::GridT.find_R2st[iat][offset];
					for(int ib=0; ib<GlobalC::pw.bxyz; ++ib					)
					{
						if(cal_flag[ib][ia1])
						{
							dgemv_(&trans, &block_size[ia1], &block_size[ia1], &alpha,
									&DMR[DM_start], &block_size[ia1], 
									&psi[ib][idx1], &inc,  
									&beta, &psi_DMR.ptr_2D[ib][idx1], &inc);
						}
					}
				}
			}

			// get (j,beta,R2)
			int start;
			switch(job)
			{
				case 1:
					start=ia1+1;
					break;
				case 2:
					start=0;
					break;
				default:
					ModuleBase::WARNING_QUIT("psi_dmr","job can only be 1 or 2");
			}

			for (int ia2=start; ia2<na_grid; ia2++)
			{
				if(all_out_of_range[ia2]) continue;

				//---------------------------------------------
				// check if we need to calculate the big cell.
				//---------------------------------------------
				bool same_flag = false;
				for(int ib=0; ib<gt.bxyz; ++ib)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						same_flag = true;
						break;
					}
				}

				if(!same_flag) continue;

				const int bcell2 = gt.bcell_start[grid_index] + ia2;
				const int T2 = GlobalC::ucell.iat2it[ gt.which_atom[bcell2]];
				const int iat2 = gt.which_atom[bcell2];

				Atom *atom2 = &GlobalC::ucell.atoms[T2];

				//---------------
				// get cell R2.
				//---------------
				const int id2 = gt.which_unitcell[bcell2];
				const int R2x = gt.ucell_index2x[id2];
				const int R2y = gt.ucell_index2y[id2];
				const int R2z = gt.ucell_index2z[id2];

				//------------------------------------------------
				// calculate the 'offset': R2 position relative
				// to R1 atom.
				//------------------------------------------------
				const int dRx = R1x - R2x;
				const int dRy = R1y - R2y;
				const int dRz = R1z - R2z;

				const int index = gt.cal_RindexAtom(dRx, dRy, dRz, iat2);
				int offset = -1;

				int* find_start = gt.find_R2[iat];
				int* findend = gt.find_R2[iat] + gt.nad[iat];

				// the nad should be a large expense of time.
				for(int* find=find_start; find < findend; find++)
				{
					if( find[0] == index )
					{
						offset = find - find_start;
						break;
					}
				}

				if(offset == -1 )
				{
					ModuleBase::WARNING_QUIT("gint_k","mult_psi_DMR wrong");
				}
				assert(offset < gt.nad[iat]);

				//--------------------------------------------------------------- 
				// what I do above is to get 'offset' for atom std::pair (iat1, iat2)
				// if I want to simplify this searching for offset,
				// I should take advantage of gt.which_unitcell.
				//--------------------------------------------------------------- 

				int cal_num=0;
   				for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    			{
        			if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        			    ++cal_num;
    			}

				if(cal_num>GlobalC::pw.bxyz/4)
				{
					const int idx1=block_index[ia1];
			        const int idx2=block_index[ia2];
    				const int DM_start = GlobalC::GridT.nlocstartg[iat]+ GlobalC::GridT.find_R2st[iat][offset];
    				dgemm_(&trans, &trans, &block_size[ia2], &GlobalC::pw.bxyz, &block_size[ia1], &alpha1,
    					&DMR[DM_start], &block_size[ia2], 
    					&psi[0][idx1], &LD_pool,
    					&beta, &psi_DMR.ptr_2D[0][idx2], &LD_pool);
				}
				else if(cal_num>0)
				{
					const int idx1=block_index[ia1];
					const int idx2=block_index[ia2];
    				const int DM_start = GlobalC::GridT.nlocstartg[iat]+ GlobalC::GridT.find_R2st[iat][offset];
					
    				for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    				{
        				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        				{
            				dgemv_(&trans, &block_size[ia2], &block_size[ia1], &alpha1,
            					&DMR[DM_start], &block_size[ia2], 
            					&psi[ib][idx1], &inc,
            					&beta, &psi_DMR.ptr_2D[ib][idx2], &inc);
        				}
    				}
				} // cal_num
			}// ia2
		}//ia1

		delete[] all_out_of_range;
		return psi_DMR;

	}
}