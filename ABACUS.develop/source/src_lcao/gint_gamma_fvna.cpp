#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"

#include "global_fp.h" // mohan add 2021-01-30

void Gint_Gamma::cal_force_vna(const double* vlocal_in, const Grid_Technique &gt, LCAO_Matrix &lm)
{
	timer::tick("Gint_Gamma","cal_force_vna",'I');
	this->vlocal = vlocal_in;
	this->save_atoms_on_grid(gt);
	this->gamma_force_vna(gt, lm);
	timer::tick("Gint_Gamma","cal_force_vna",'I');
}


//-------------------------------------------------------------
// Because this part calculat <dphi/dR | Vna | phi>
// so there are many reasons make this part very expensive.
// (1) Dense grid. (dx*2,dy*2,dz*2) 8 times at least
// (2) All the matrix are needed, it's non-symmetric. 2 times at least
// (3) Need to calculate three components of force: Fx,Fy,Fz
// However, 3 times may be
// 8*2*3=48 (20)
// we must know,
// once we have the lattice structure and the atom positions,
// and of course the local orbitals.
// we can do the neutral potential force right away,
// why should we bother so many processors to do at the same
// time after the charge density is converged.
// we can begin to do it in another set of processors!!
// 
// The same principle applys to the <phi|Vna|phi>
//-------------------------------------------------------------
void Gint_Gamma::gamma_force_vna(const Grid_Technique &gt, LCAO_Matrix &lm)
{
    TITLE("Grid_Integral","gamma_force_vna");
    timer::tick("Gint_Gamma","gamma_force_vna",'J');
	// gt.lgd: local grid dimension (sub-FFT-mesh).
    double** DGridV = new double*[gt.lgd];
    double** DGridV_s;
    if(STRESS) DGridV_s = new double*[gt.lgd];
    for (int i=0; i<gt.lgd; ++i)
    {
        DGridV[i] = new double[gt.lgd*3];
       	ZEROS(DGridV[i], 3*gt.lgd);
        if(STRESS)
        {
            DGridV_s[i] = new double[gt.lgd*6];
            ZEROS(DGridV_s[i], 6*gt.lgd);
        }
    }
    Memory::record("Gint_Gamma","DGridV",3*gt.lgd*gt.lgd,"double");
    if(STRESS) Memory::record("Gint_Gamma","DGridV_s",6*gt.lgd*gt.lgd,"double");
    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;
	if( delta_r != ORB.Vna[0].dr )
	{
		cout << " delta_r = " << delta_r << endl;
		cout << " dr = " << ORB.Vna[0].dr << endl;
		WARNING_QUIT("Gint_Gamma::gamma_force_vna","delta_r != ORB.Vna[0].dr");
	}

	const int bxyz = gt.bxyz;
	const int bx=gt.bx;
	const int by=gt.by;
	const int bz=gt.bz;

    double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
    double** distance; // distance between atom and grid: [bxyz, maxsize]
    double*** psir_ylm;
	double*** dphi;
    bool** cal_flag;
    if(max_size!=0)
    {
        dr = new double**[bxyz];
        distance = new double*[bxyz];
        psir_ylm = new double**[bxyz];
		dphi = new double**[bxyz];
        cal_flag = new bool*[bxyz];
        for(int i=0; i<bxyz; i++)
        {
            dr[i] = new double*[max_size];
            distance[i] = new double[max_size];
            psir_ylm[i] = new double*[max_size];
			dphi[i] = new double*[max_size];
            cal_flag[i] = new bool[max_size];

            ZEROS(distance[i], max_size);
            ZEROS(cal_flag[i], max_size);

            for(int j=0; j<max_size; j++)
            {
                dr[i][j] = new double[3];
                psir_ylm[i][j] = new double[ucell.nwmax];
                dphi[i][j] = new double[ucell.nwmax*3]; //3 includes x,y,z
                ZEROS(dr[i][j],3);
                ZEROS(psir_ylm[i][j],ucell.nwmax);
				ZEROS(dphi[i][j],ucell.nwmax*3);
            }
        }
    }

    int nnnmax=0;
    for(int T=0; T<ucell.ntype; T++)
    {
        nnnmax = max(nnnmax, nnn[T]);
    }

    //array to store spherical harmonics and its derivatives
	assert(nnnmax<400);
	// Peize Lin change rly, grly 2016-08-26
    vector<double> rly;
    vector<vector<double>> grly;

    double mt[3]={0,0,0};
    double *vldr3 = new double[bxyz];
    ZEROS(vldr3, bxyz);

	double* vna3d = new double[bxyz];
	ZEROS(vna3d, bxyz);


	int* nww = new int[max_size];
	int* iw0_all = new int[max_size];
	ZEROS(nww,max_size);
	ZEROS(iw0_all,max_size);

//	ofstream ofs_x("vna_x.dat");
//	ofstream ofs_y("vna_y.dat");
//	ofstream ofs_z("vna_z.dat");

	const int nbx = gt.nbx;
	const int nby = gt.nby;
	const int nbz_start = gt.nbzp_start;
	const int nbz = gt.nbzp;

	
	// needed in the inner circle.
	double vv;
	double* psi1;
	double* dphi2;
	double* dpvp;
        double* dpvpr;
	double* dphi3;
	double* phi2_end;
	int iw1_lo;

    for (int i=0; i<nbx; i++)
    {
        for (int j=0; j<nby; j++)
        {
            for (int k=nbz_start; k<nbz_start+nbz; k++)
            {
                this->grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
                const int size = gt.how_many_atoms[ this->grid_index ];
				if(size==0)continue; // tmp by mohan
					
				ZEROS(vna3d, bxyz);//vna

                for (int id=0; id<size; id++)
                {
                    const int mcell_index = gt.bcell_start[grid_index] + id;
                    const int imcell = gt.which_bigcell[mcell_index];
                    int iat = gt.which_atom[mcell_index];
                    const int it = ucell.iat2it[ iat ];
                    const int ia = ucell.iat2ia[ iat ];
                    Atom *atom = &ucell.atoms[it];

					mt[0] = gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0];
					mt[1] = gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1];
					mt[2] = gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2];

					for(int ib=0; ib<bxyz; ib++)
                    {
                        dr[ib][id][0] = gt.meshcell_pos[ib][0] + mt[0];
                        dr[ib][id][1] = gt.meshcell_pos[ib][1] + mt[1];
                        dr[ib][id][2] = gt.meshcell_pos[ib][2] + mt[2];

                        distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
						+ dr[ib][id][1]*dr[ib][id][1] 
						+ dr[ib][id][2]*dr[ib][id][2]);
                        
						if(distance[ib][id] > ORB.Vna[it].rcut && 
							distance[ib][id] > ORB.Phi[it].getRcut()) 
						{
							cal_flag[ib][id]=false;
							continue;
						}
						
                    	Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[ib][id][0], 
						dr[ib][id][1], dr[ib][id][2], rly, grly);

                    	// 1E-7 is necessary in case of R is just on one grid
                    	const double position = distance[ib][id] / delta_r;
                    	this->iq[id] = static_cast<int>(position);
                    	this->x0[id] = position - static_cast<double>(iq[id]);
                    	this->x1[id] = 1.0 - x0[id];
                    	this->x2[id] = 2.0 - x0[id];
                    	this->x3[id] = 3.0 - x0[id];
                    	this->x12[id] = x1[id]*x2[id] / 6.0;
						this->x03[id] = x0[id]*x3[id] / 2.0;
						

						// vna						
						if(distance[ib][id] <= ORB.Vna[it].rcut)
						{
							double ccc = (x12[id]*(ORB.Vna[it].vna_u[iq[id]]*x3[id]
							 +ORB.Vna[it].vna_u[iq[id]+3]*x0[id])
							+ x03[id]*(ORB.Vna[it].vna_u[iq[id]+1]*x2[id]
							-ORB.Vna[it].vna_u[iq[id]+2]*x1[id])); 

							vna3d[ib] += ccc;
							
//							cout << " size=" << size << " distance=" << distance[ib][id] << " ccc=" << ccc << endl;

//							vna3d[ib] += (x12[id]*(ORB.Vna[it].vna_u[iq[id]]*x3[id]
//										+ORB.Vna[it].vna_u[iq[id]+3]*x0[id])
//									+ x03[id]*(ORB.Vna[it].vna_u[iq[id]+1]*x2[id]
//										-ORB.Vna[it].vna_u[iq[id]+2]*x1[id]));/// sqrt(FOUR_PI);
						}


						if(distance[ib][id] <= ORB.Phi[it].getRcut())
                        {
                            cal_flag[ib][id]=true;
                        }
                        else
                        {
                            cal_flag[ib][id]=false;
                            continue;
                        }


						
						double tmp, dtmp;
						for (int iw=0; iw< atom->nw; iw++)
						{
							if ( atom->iw2_new[iw] )
							{
								pointer = &ORB.Phi[it].PhiLN(
										atom->iw2l[iw],
										atom->iw2n[iw]);

								if ( iq[id] >= pointer->nr_uniform-4)
								{
									tmp = dtmp = 0.0;
								}
								else
								{
									tmp = x12[id]*(pointer->psi_uniform[iq[id]]*x3[id]
											+pointer->psi_uniform[iq[id]+3]*x0[id])
										+ x03[id]*(pointer->psi_uniform[iq[id]+1]*x2[id]
												-pointer->psi_uniform[iq[id]+2]*x1[id]);

									dtmp = x12[id]*(pointer->dpsi_uniform[iq[id]]*x3[id]
											+pointer->dpsi_uniform[iq[id]+3]*x0[id])
										+ x03[id]*(pointer->dpsi_uniform[iq[id]+1]*x2[id]
												-pointer->dpsi_uniform[iq[id]+2]*x1[id]);
								}
							}//new l is used.

							int ll = atom->iw2l[iw];
							int idx_lm = atom->iw2_ylm[iw];

							//special case for distance[id] -> 0
							//Problems Remained
							//You have to add this two lines
							double rr = distance[ib][id];
							if (rr < 1e-9)
							{
								if (ll == 0)
								{
									psir_ylm[ib][id][iw] = tmp * rly[idx_lm];
									dphi[ib][id][3*iw] = dphi[ib][id][3*iw+1] = dphi[ib][id][3*iw+2] = 0.0;
								}
								else
								{
									pointer = &ORB.Phi[it].
										PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

									double Zty = pointer->zty;
									psir_ylm[ib][id][iw] = Zty * rly[idx_lm];

									const int iww = 3*iw; 
									dphi[ib][id][iww] = Zty * grly[idx_lm][0];
									dphi[ib][id][iww+1] = Zty * grly[idx_lm][1];
									dphi[ib][id][iww+2] = Zty * grly[idx_lm][2];
								}
							}
							else
							{
								double rl;
								if(ll==0)
								{
									rl = 1.0;
								}
								else if(ll==1)
								{
									rl = rr;
								}
								else 
								{
									rl = pow(rr, ll);
								}

								psir_ylm[ib][id][iw] = tmp * rly[idx_lm] / rl;

								double tmpdphi_rly = (dtmp  - tmp * ll / rr) / rl * rly[idx_lm] / rr;
								double tmprl = tmp/rl;

								dphi[ib][id][iw*3] = tmpdphi_rly * dr[ib][id][0]  + tmprl * grly[idx_lm][0];
								dphi[ib][id][iw*3+1] = tmpdphi_rly * dr[ib][id][1]  + tmprl * grly[idx_lm][1];
								dphi[ib][id][iw*3+2] = tmpdphi_rly * dr[ib][id][2]  + tmprl * grly[idx_lm][2];
							}
						}
					}// ib
				}//!id //finish loop of calc pre-info for each adjacent atom

				int bindex = 0;
				// z is the fastest,
				for(int ii=0; ii<bx; ii++)
				{
					for(int jj=0; jj<by; jj++)
					{
						for(int kk=0; kk<bz; kk++)
						{
							// 0.074 is compensation for Si.
							//vldr3[bindex] = ( vna3d[bindex]+0.074)  * this->vfactor;
							vldr3[bindex] = vna3d[bindex]  * this->vfactor;
							++bindex;
						}
					}
				}



				for(int ia1=0; ia1<size; ++ia1)
				{
					const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
					const int iat1 = gt.which_atom[mcell_index1];
					const int T1 = ucell.iat2it[iat1];
					const int I1 = ucell.iat2ia[iat1];
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
					iw0_all[ia1] = gt.trace_lo[start1];
					nww[ia1] = ucell.atoms[T1].nw;
				}

				// The first method
				// make pair among adjacent atoms of grid with number size
				//(ia1, ia2)
				for (int ia1=0; ia1<size; ++ia1)
				{
					for (int ia2=0; ia2<size; ++ia2)
					{
						const int iw2_all = 3*iw0_all[ia2];
						const int iw2_s = 6*iw0_all[ia2];
						const int nww2 = 3*nww[ia2];
						for(int ib=0; ib<bxyz; ++ib)
						{
							const double vnow = vldr3[ib];
							if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
							{
								psi1 = &psir_ylm[ib][ia1][0];
								dphi2 = dphi[ib][ia2];
								iw1_lo = iw0_all[ia1];
								// how many orbitals in this type: SZ or DZP or TZP...
								for (int iw=0; iw< nww[ia1]; ++iw, ++psi1)
								{
									vv = vnow * psi1[0];
									dpvp = &DGridV[iw1_lo+iw][iw2_all];
									if(STRESS) dpvpr = &DGridV_s[iw1_lo+iw][iw2_s];
									dphi3 = dphi2;
									phi2_end = dphi3 + nww2;
									while(dphi3<phi2_end)
									{

//										dpvp[0] -= dphi3[0] * vv;
//										++dpvp;
//										++dphi3;
//										dpvp[0] -= dphi3[0] * vv;
//										++dpvp;
//										++dphi3;	
//										dpvp[0] -= dphi3[0] * vv;
//										++dpvp;
//										++dphi3;

										dpvp[0] -= dphi3[0] * vv;
										dpvp[1] -= dphi3[1] * vv;
										dpvp[2] -= dphi3[2] * vv;
										dpvp+=3;
										if(STRESS)
										{
											dpvpr[0] -= dphi3[0] * vv * dr[ib][ia2][0];
											dpvpr[1] -= dphi3[0] * vv * dr[ib][ia2][1];
											dpvpr[2] -= dphi3[0] * vv * dr[ib][ia2][2];
											dpvpr[3] -= dphi3[1] * vv * dr[ib][ia2][1];
											dpvpr[4] -= dphi3[1] * vv * dr[ib][ia2][2];
											dpvpr[5] -= dphi3[2] * vv * dr[ib][ia2][2];
											dpvpr+=6;
										}
										dphi3+=3;
									}
								}//iw1
							}// cal_flag
						}//ib
					}//!ia
				}// ia1

				// The second method
				// make pair among adjacent atoms of grid with number size
				//(ia1, ia2)
				/*
				for (int ia1=0; ia1<size; ++ia1)
				{
					for (int ia2=0; ia2<size; ++ia2)
					{
						const int iw2_all = 3*iw0_all[ia2];
						int iw1_lo = iw0_all[ia1];
						// how many orbitals in this type: SZ or DZP or TZP...
						for (int iw=0; iw< nww[ia1]; ++iw)
						{
							double* dpvp = &DGridV[iw1_lo][iw2_all];
							for(int iw2=0; iw2<nww[ia2]; ++iw2)
							{
								const int iw22 = iw2*3;
								double vvx = 0.0;
								double vvy = 0.0;
								double vvz = 0.0;
								for(int ib=0; ib<bxyz; ++ib)
								{
									if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
									{
										vvx += psir_ylm[ib][ia1][iw] * vldr3[ib] * dphi[ib][ia2][iw22];
										vvy += psir_ylm[ib][ia1][iw] * vldr3[ib] * dphi[ib][ia2][iw22+1];
										vvz += psir_ylm[ib][ia1][iw] * vldr3[ib] * dphi[ib][ia2][iw22+2];
									}
								}
								dpvp[0] -= vvx;
								++dpvp;
								dpvp[0] -= vvy;
								++dpvp;
								dpvp[0] -= vvz;
								++dpvp;
							}
							iw1_lo++;
						}//iw1
					}//!ia
				}// ia1
				*/

				// The third method
				// make pair among adjacent atoms of grid with number size
				//(ia1, ia2)
				/*
				for (int ia1=0; ia1<size; ++ia1)
				{
					for(int ib=0; ib<bxyz; ++ib)
					{
						if(cal_flag[ib][ia1])
						{
							psi1 = &psir_ylm[ib][ia1][0];
							iw1_lo = iw0_all[ia1];
							// how many orbitals in this type: SZ or DZP or TZP...
							for (int iw=0; iw< nww[ia1]; ++iw, ++psi1)
							{
								vv = vldr3[ib] * psi1[0];
								for (int ia2=0; ia2<size; ++ia2)
								{
									if(cal_flag[ib][ia2])
									{
										dphi2 = dphi[ib][ia2];
										dpvp = &DGridV[iw1_lo][3*iw0_all[ia2]];
										phi2_end = dphi2 + 3*nww[ia2];
										while(dphi2<phi2_end)
										{
											// this is not efficient
//											dpvp[0] -= dphi2[0] * vv;
//											++dpvp;
//											++dphi2;
//											dpvp[0] -= dphi2[0] * vv;
//											++dpvp;
//											++dphi2;	
//											dpvp[0] -= dphi2[0] * vv;
//											++dpvp;
//											++dphi2;
											// this is efficient!
											dpvp[0] -= dphi2[0] * vv;
											dpvp[1] -= dphi2[1] * vv;
											dpvp[2] -= dphi2[2] * vv;
											dpvp+=3;
											dphi2+=3;
										}
									}
								}
								++iw1_lo;
							}//iw
						}//cal_flag
					}//ib
				}// ia1
				*/


			}// k
		}// j
	}// i

	delete[] vna3d;
	delete[] nww;
	delete[] iw0_all;

	
    if(max_size!=0)
    {
        for(int i=0; i<bxyz; i++)
        {
            for(int j=0; j<max_size; j++)
            {
                delete[] dr[i][j];
                delete[] psir_ylm[i][j];
				delete[] dphi[i][j];
            }
            delete[] dr[i];
            delete[] distance[i];
            delete[] psir_ylm[i];
            delete[] cal_flag[i];
			delete[] dphi[i];
        }
        delete[] dr;
        delete[] distance;
        delete[] psir_ylm;
        delete[] cal_flag;
		delete[] dphi;
    }
	delete[] vldr3;

#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
    timer::tick("Gint_Gamma","gamma_force_vna",'J');

    timer::tick("Gint_Gamma","distri_fvna",'J');
	//-------------------------------------------
	// distribution of the Hamiltonian matrix.
	//-------------------------------------------
    double* tmp = new double[NLOCAL*3];
    double* tmpr;
    if(STRESS) tmpr = new double[NLOCAL*6];
    for (int i=0; i<NLOCAL; i++)
    { 
		ZEROS(tmp, 3*NLOCAL);   
		if(STRESS) ZEROS(tmpr, 6*NLOCAL);
		const int mu = gt.trace_lo[i];
		// mohan fix bug 2010-09-05
		// lack mu>=0 and nu>=0 in previous version.
		if(mu >=0)
		{
			for (int j=0; j<NLOCAL; j++)
			{
				const int nu = gt.trace_lo[j];
				if(nu>=0)
				{
					tmp[3*j]   = DGridV[mu][3*nu];
					tmp[3*j+1] = DGridV[mu][3*nu+1];
					tmp[3*j+2] = DGridV[mu][3*nu+2];
					if(STRESS) 
					{
						tmpr[6*j]   = DGridV_s[mu][6*nu];
						tmpr[6*j+1] = DGridV_s[mu][6*nu+1];
						tmpr[6*j+2] = DGridV_s[mu][6*nu+2];
						tmpr[6*j+3] = DGridV_s[mu][6*nu+3];
						tmpr[6*j+4] = DGridV_s[mu][6*nu+4];
						tmpr[6*j+5] = DGridV_s[mu][6*nu+5];
					}
				}
			}
		}

		Parallel_Reduce::reduce_double_pool( tmp, 3*NLOCAL );
		if(STRESS) Parallel_Reduce::reduce_double_pool( tmpr, 6*NLOCAL );
		const int i2d = ParaO.trace_loc_row[i];
		if(i2d<0) continue;
		for (int j=0; j<NLOCAL; j++)
		{
			const int j2d = ParaO.trace_loc_col[j];
			if(j2d<0) continue;
			const int index = i2d * ParaO.ncol + j2d; 
		   	LM.DHloc_fixed_x[index] += tmp[3*j];
			LM.DHloc_fixed_y[index] += tmp[3*j+1];
			LM.DHloc_fixed_z[index] += tmp[3*j+2];
			if(STRESS)
			{
				LM.DHloc_fixed_11[index] += tmpr[6*j];
				LM.DHloc_fixed_12[index] += tmpr[6*j+1];
				LM.DHloc_fixed_13[index] += tmpr[6*j+2];
				LM.DHloc_fixed_22[index] += tmpr[6*j+3];
				LM.DHloc_fixed_23[index] += tmpr[6*j+4];
				LM.DHloc_fixed_33[index] += tmpr[6*j+5];
			}
//            lm.set_force (i,j,tmp[j], tmp[j+1], tmp[j+2],'N');
		}
	}
	delete[] tmp;
	if(STRESS) delete[] tmpr;
	//delete DGridV_x,y,z
	for (int i=0; i<gt.lgd; ++i)
	{
		delete[] DGridV[i];
		if(STRESS) delete[] DGridV_s[i];
	}
	delete[] DGridV;
	if(STRESS) delete[] DGridV_s;

    timer::tick("Gint_Gamma","distri_fvna",'J');
    return;
}

