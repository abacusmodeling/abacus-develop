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
	const double *vl)
{
	ModuleBase::TITLE("Gint_k","cal_force_k");
	ModuleBase::timer::tick("Gint_k","cal_force_k");

	int nnrg = GlobalC::GridT.nnrg;

	if(GlobalV::OUT_LEVEL != "m") GlobalV::ofs_running << " LNNR.nnrg in cal_force_k = " << nnrg << std::endl;
	assert(nnrg>=0);

	// just because to make thea arrys meaningful.
	if(nnrg == 0)
	{
		nnrg = 1;
	}

	// to store < phi | vlocal | dphi>
	double* pvdpx;
	double* pvdpy;
	double* pvdpz;
	double* pvdp11;
	double* pvdp22;
	double* pvdp33;
	double* pvdp12;
	double* pvdp13;
	double* pvdp23;
	if(isforce)
	{
		pvdpx = new double[nnrg];
		pvdpy = new double[nnrg];
		pvdpz = new double[nnrg];
		ModuleBase::GlobalFunc::ZEROS(pvdpx, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdpy, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdpz, nnrg);
	}
	if(isstress)
	{
		pvdp11 = new double[nnrg];
		pvdp22 = new double[nnrg];
		pvdp33 = new double[nnrg];
		pvdp12 = new double[nnrg];
		pvdp13 = new double[nnrg];
		pvdp23 = new double[nnrg];
		ModuleBase::GlobalFunc::ZEROS(pvdp11, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdp22, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdp33, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdp12, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdp13, nnrg);
		ModuleBase::GlobalFunc::ZEROS(pvdp23, nnrg);
	}

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

					if(isforce)
					{
						this->cal_meshball_force(grid_index, na_grid,
							block_index,
							cal_flag, psir_vlbr3.ptr_2D,
							dpsir_ylm_x.ptr_2D,
							dpsir_ylm_y.ptr_2D,
							dpsir_ylm_z.ptr_2D,
							pvdpx, pvdpy, pvdpz,
							GlobalC::GridT);
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
						this->cal_meshball_stress(grid_index, na_grid,
							block_index,
							cal_flag, psir_vlbr3.ptr_2D,
                            dpsir_ylm_xx.ptr_2D, 
                            dpsir_ylm_xy.ptr_2D, 
                            dpsir_ylm_xz.ptr_2D,
                            dpsir_ylm_yy.ptr_2D, 
                            dpsir_ylm_yz.ptr_2D, 
                            dpsir_ylm_zz.ptr_2D,
							pvdp11, pvdp22, pvdp33, pvdp12, pvdp13, pvdp23,
							GlobalC::GridT);
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

	//---------------------------------------
	// Folding R here
	//---------------------------------------

	//this->LM->DHloc_fixedR_x
	this->folding_force(isforce, isstress, fvl_dphi, svl_dphi, pvdpx, pvdpy, pvdpz,
			pvdp11, pvdp22, pvdp33, pvdp12, pvdp13, pvdp23);

	if(isforce)
	{
		delete[] pvdpx;
		delete[] pvdpy;
		delete[] pvdpz;
	}
	if(isstress)
	{
		delete[] pvdp11;
		delete[] pvdp22;
		delete[] pvdp33;
		delete[] pvdp12;
		delete[] pvdp13;
		delete[] pvdp23;
	}

	ModuleBase::timer::tick("Gint_k","cal_force_k");
	return;
}

#include "../module_base/mathzone.h"
void Gint_k::cal_meshball_stress(
	const int &grid_index, 
	const int &size,
	const int*const block_index,
	bool** cal_flag, 
	double** psir_vlbr3,
    double** dpsir_xx,
    double** dpsir_xy,
    double** dpsir_xz,
    double** dpsir_yy,
    double** dpsir_yz,
    double** dpsir_zz,
	double* pvdp11, 
	double* pvdp22, 
	double* pvdp33, 
	double* pvdp12, 
	double* pvdp13, 
	double* pvdp23, 
	const Grid_Technique &gt)
{
	ModuleBase::timer::tick("Gint_k","evaluate_vl_stress");
	double *psi2;
    double *pvp11, *pvp22, *pvp33, *pvp12, *pvp13, *pvp23;
	int iwi, iww;
	double *psixx, *psixy, *psixz, *psiyy, *psiyz, *psizz;

	bool *all_out_of_range = new bool[size];
	for(int ia=0; ia<size; ++ia) //number of atoms
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

	double* dmR = this->DM_R[GlobalV::CURRENT_SPIN];

	double* dmR2;
	for (int ia1=0; ia1<size; ++ia1)
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
        const int DM_start = gt.nlocstartg[iat];

        // get (j,beta,R2)
        for (int ia2=0; ia2<size; ++ia2)
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
				//	std::cout << " same flag is = " << ib << std::endl;
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
				ModuleBase::WARNING_QUIT("gint_k","evaluate_vl_force wrong");
			}
			assert(offset < gt.nad[iat]);

			//--------------------------------------------------------------- 
			// what I do above is to get 'offset' for atom std::pair (iat1, iat2)
			// if I want to simplify this searching for offset,
			// I should take advantage of gt.which_unitcell.
			//--------------------------------------------------------------- 

			const int iatw = DM_start + gt.find_R2st[iat][offset];
			
			for(int ib=0; ib<gt.bxyz; ++ib)
			{
				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
				{
					psi2 = &psir_vlbr3[ib][block_index[ia2]];

					psixx = &dpsir_xx[ib][block_index[ia1]];
					psixy = &dpsir_xy[ib][block_index[ia1]];
					psixz = &dpsir_xz[ib][block_index[ia1]];
					psiyy = &dpsir_yy[ib][block_index[ia1]];
					psiyz = &dpsir_yz[ib][block_index[ia1]];
					psizz = &dpsir_zz[ib][block_index[ia1]];
						
					//------------------------------------
					// circle for wave functions of atom 1.
					//------------------------------------
					iwi = 0;

					for (int iw1=0; iw1 < atom1->nw; iw1++)
					{

						iww = iatw + iwi;// -1 because ++iww from below.
						dmR2 = &dmR[iww]; //mohan add 2012-01-05
						iwi += atom2->nw;

						pvp11 = &pvdp11[iww]; //zhengdy add 2017/3/28
						pvp22 = &pvdp22[iww];
						pvp33 = &pvdp33[iww];
						pvp12 = &pvdp12[iww]; 
						pvp13 = &pvdp13[iww];
						pvp23 = &pvdp23[iww];
						//------------------------------------
						// circle for wave functions of atom 2.
						//------------------------------------
						for(int iw2=0; iw2 < atom2->nw; iw2++)
						{

							pvp11[iw2] += dmR2[iw2] * psixx[iw1] * psi2[iw2];
							pvp22[iw2] += dmR2[iw2] * psiyy[iw1] * psi2[iw2];
							pvp33[iw2] += dmR2[iw2] * psizz[iw1] * psi2[iw2];
							pvp12[iw2] += dmR2[iw2] * psixy[iw1] * psi2[iw2];
							pvp13[iw2] += dmR2[iw2] * psixz[iw1] * psi2[iw2];
							pvp23[iw2] += dmR2[iw2] * psiyz[iw1] * psi2[iw2];
						}
					}// iw
				}//end flag
			}//end ib
        }// ia2
	}//ia1

	delete[] all_out_of_range;
	ModuleBase::timer::tick("Gint_k","evaluate_vl_stress");
	return;
}


void Gint_k::cal_meshball_force(const int &grid_index, const int &size,
	const int*const block_index,
	bool** cal_flag, double** psir_vlbr3, 
	double** dphi_x, double** dphi_y, double** dphi_z,
	double* pvdpx, double* pvdpy, double* pvdpz,
	const Grid_Technique &gt)
{                       
	ModuleBase::timer::tick("Gint_k","evaluate_vl_force");
	double *psi2;
	double *pvp1, *pvp2, *pvp3;//extra pointer
	int iwi, iww;
	double *psix, *psiy, *psiz;

	bool *all_out_of_range = new bool[size];
	for(int ia=0; ia<size; ++ia) //number of atoms
	{
		all_out_of_range[ia] = true;
		for(int ib=0; ib<gt.bxyz; ++ib) //number of small box in big box
		{
			if(cal_flag[ib][ia])
			{
				all_out_of_range[ia] = false;
			}
		}
	}

	double* dmR = this->DM_R[GlobalV::CURRENT_SPIN];
	double* dmR2;
	for (int ia1=0; ia1<size; ia1++)
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
		const int DM_start = gt.nlocstartg[iat];

		// get (j,beta,R2)
		for (int ia2=0; ia2<size; ia2++)
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
				ModuleBase::WARNING_QUIT("gint_k","evaluate_vl_force wrong");
			}
			assert(offset < gt.nad[iat]);

			//--------------------------------------------------------------- 
			// what I do above is to get 'offset' for atom std::pair (iat1, iat2)
			// if I want to simplify this searching for offset,
			// I should take advantage of gt.which_unitcell.
			//--------------------------------------------------------------- 

			const int iatw = DM_start + gt.find_R2st[iat][offset];

			for(int ib=0; ib<gt.bxyz; ++ib)
			{
				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
				{
					psi2 = &psir_vlbr3[ib][block_index[ia2]];

					psix = &dphi_x[ib][block_index[ia1]];
					psiy = &dphi_y[ib][block_index[ia1]];
					psiz = &dphi_z[ib][block_index[ia1]];
					//------------------------------------
					// circle for wave functions of atom 1.
					//------------------------------------
					iwi = 0;

					for (int iw1=0; iw1 < atom1->nw; iw1++)
					{
						iww = iatw + iwi;// -1 because ++iww from below.
						dmR2 = &dmR[iww]; //mohan add 2012-01-05
						iwi += atom2->nw;

						//--------------------------------------
						// store phi_i(r) vlocal(r) * dphi_j(r)
						//--------------------------------------
						pvp1 = &pvdpx[iww]; //mohan add 2012-1-6
						pvp2 = &pvdpy[iww];
						pvp3 = &pvdpz[iww];

						//------------------------------------
						// circle for wave functions of atom 2.
						//------------------------------------
						for(int iw2=0; iw2 < atom2->nw; iw2++)
						{
							// the main difference to calculate
							// the force is that the whole 
							// matrix should be calculated!

							// DM(R) * psi1(r) * v(r) * psi2_R(r)
							pvp1[iw2] += dmR2[iw2] * psix[iw1] * psi2[iw2];
							pvp2[iw2] += dmR2[iw2] * psiy[iw1] * psi2[iw2];
							pvp3[iw2] += dmR2[iw2] * psiz[iw1] * psi2[iw2];
							//density matrix
						}
					}// iw
				}//end flag
			}//end ib
		}// ia2
	}//ia1

	delete[] all_out_of_range;
	ModuleBase::timer::tick("Gint_k","evaluate_vl_force");
	return;
}