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

	Array_Pool<double> cal_psir_ylm(
		const int na_grid, 					// number of atoms on this grid 
		const int LD_pool,
		const int grid_index, 				// 1d index of FFT index (i,j,k) 
		const double delta_r, 				// delta_r of the uniform FFT grid
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag) 	// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	{
		Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);
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
				double *p=&psir_ylm.ptr_2D[ib][block_index[id]];
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
		return psir_ylm;
	}

}