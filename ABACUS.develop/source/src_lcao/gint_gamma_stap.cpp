//sun zhiyuan add at 2011-08-10, for calculating A*P for bfield//
#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "bfield.h"
void Gint_Gamma::cal_S_T_AP(char type, const Grid_Technique &gt)
{
	this->save_atoms_on_grid(gt);
    assert(bfid.Ready_A_of_Atom==1);
	this->gamma_S_T_AP(type, gt);
}

// this subroutine lies in the heart of LCAO algorithms.
// so it should be done very efficiently, very carefully.
// I might repeat again to emphasize this: need to optimize
// this code very efficiently, very carefully.
void Gint_Gamma::gamma_S_T_AP(char type, const Grid_Technique &gt)
{
	TITLE("Grid_Integral","gamma_S_T_AP");
	timer::tick("Gint_Gamma","gamma_S_T_AP");

	cout << " calculate " << type << endl;

	//get cartesian lattice vectors, in bohr//
	double latvec1[3],latvec2[3],latvec3[3];
	latvec1[0]=ucell.latvec.e11*ucell.lat0;
	latvec1[1]=ucell.latvec.e12*ucell.lat0;
	latvec1[2]=ucell.latvec.e13*ucell.lat0;
	latvec2[0]=ucell.latvec.e21*ucell.lat0;
	latvec2[1]=ucell.latvec.e22*ucell.lat0;
	latvec2[2]=ucell.latvec.e23*ucell.lat0;
	latvec3[0]=ucell.latvec.e31*ucell.lat0;
	latvec3[1]=ucell.latvec.e32*ucell.lat0;
	latvec3[2]=ucell.latvec.e33*ucell.lat0;

	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	const double delta_r = ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;
	double mt[3]={0,0,0};
	double A_atom1[3]={0,0,0};
    double A_atom2[3]={0,0,0};
    double A2_ms_A1[3]={0,0,0}; //used to store A2-A1//
	//-----------------------------------------------------------
	//temporary matrixes of S, T, A*P and A*A. lgd is the number of related orbitals//  
	complex<double> ** GridS;  //for S matrix//
	complex<double> ** Grid_T;  //for T matrix//
	complex<double> ** GridAP; //for A*P matrix//
	complex<double> ** GridAA; //for A*A matrix//
    bool perform_gint = true;
	const int lgd_now = gt.lgd;
	if(lgd_now > 0)
	{
		GridS = new complex<double>* [lgd_now];
		Grid_T = new complex<double>* [lgd_now];
		GridAP = new complex<double>* [lgd_now];
		GridAA = new complex<double>* [lgd_now];
		for (int i=0; i<lgd_now; i++)
		{
			GridS[i] = new complex<double> [lgd_now];
			Grid_T[i] = new complex<double> [lgd_now];
			GridAP[i] = new complex<double> [lgd_now];
			GridAA[i] = new complex<double> [lgd_now];
			ZEROS(GridS[i], lgd_now);
			ZEROS(Grid_T[i], lgd_now);
			ZEROS(GridAP[i], lgd_now);
			ZEROS(GridAA[i], lgd_now);
		}
		Memory::record("Gint_Gamma","Grid_STAP",4*2*lgd_now*lgd_now,"double");
	}
	else if(lgd_now <= 0)
	{
		perform_gint = false;
	}
	//-----------------------------------------------------------
   
    //some arrays//	
	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm; //value of wave functions//
	double*** dphix; //derivatives of wave functions// 
	double*** dphiy;
	double*** dphiz;
	bool** cal_flag; //use this to decide whether the distance is larger than cutoff//
	
	if(max_size!=0)
	{
		dr = new double**[pw.bxyz];
		distance = new double*[pw.bxyz];
		psir_ylm = new double**[pw.bxyz];
		dphix = new double**[pw.bxyz];
		dphiy = new double**[pw.bxyz];
		dphiz = new double**[pw.bxyz];
		cal_flag = new bool*[pw.bxyz];
		for(int i=0; i<pw.bxyz; i++)
		{
			dr[i] = new double*[max_size];
			distance[i] = new double[max_size];
			psir_ylm[i] = new double*[max_size];
			dphix[i] = new double*[max_size];
			dphiy[i] = new double*[max_size];
			dphiz[i] = new double*[max_size];
			cal_flag[i] = new bool[max_size];

			ZEROS(distance[i], max_size);
			ZEROS(cal_flag[i], max_size);

			for(int j=0; j<max_size; j++)
			{
				dr[i][j] = new double[3];
				psir_ylm[i][j] = new double[ucell.nwmax];
				dphix[i][j] = new double[ucell.nwmax];
				dphiy[i][j] = new double[ucell.nwmax];
				dphiz[i][j] = new double[ucell.nwmax];
				ZEROS(dr[i][j],3);
				ZEROS(psir_ylm[i][j],ucell.nwmax);
				ZEROS(dphix[i][j],ucell.nwmax);
				ZEROS(dphiy[i][j],ucell.nwmax);
				ZEROS(dphiz[i][j],ucell.nwmax);
			}
		}
	}

	//array to store spherical harmonics and its derivatives
	int nnnmax=0;
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax = max(nnnmax, nnn[T]);
	}
	assert(nnnmax<400);
	// Peize Lin change rly, grly 2016-08-26
	vector<double> rly;
	vector<vector<double>> grly;

    //arrays to store the vector potential A and position r, updated at 2012-01-05//
	double **Aldr3 = new double *[pw.bxyz]; // bxyz*3//
	double **Rbxyz = new double *[pw.bxyz]; // bxyz*3//
	double *AAdr = new double [pw.bxyz];    // bxyz//
	ZEROS(AAdr, pw.bxyz);
	for(int i=0;i<pw.bxyz;i++) 
	{
		Aldr3[i]=new double[3];
		Rbxyz[i]=new double[3];
		ZEROS(Aldr3[i], 3);
		ZEROS(Rbxyz[i], 3);
	}

	const int bxyz = gt.bxyz;
	const int bx = gt.bx;
	const int by = gt.by;
	const int bz = gt.bz;

	const int nbx = gt.nbx;
	const int nby = gt.nby;
	const int nbz_start = gt.nbzp_start;
	const int nbz = gt.nbzp;

//=======================big cycle starts==============================================
	for (int i=0; i< nbx; i++)
	{
		for (int j=0; j< nby; j++)
		{
			for (int k= nbz_start; k< nbz_start+nbz; k++)
			{
				this->grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
				const int size = gt.how_many_atoms[ this->grid_index ];
				if(size==0)continue;
				
				//get pre-information about adjacent atoms//
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

					for(int ib=0; ib<pw.bxyz; ib++)
					{
						dr[ib][id][0] = gt.meshcell_pos[ib][0] + mt[0];
						dr[ib][id][1] = gt.meshcell_pos[ib][1] + mt[1];
						dr[ib][id][2] = gt.meshcell_pos[ib][2] + mt[2];

						distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
								+ dr[ib][id][1]*dr[ib][id][1] 
								+ dr[ib][id][2]*dr[ib][id][2]);

						if(distance[ib][id] <= ORB.Phi[it].getRcut())
						{
							cal_flag[ib][id]=true;
						}
						else
						{
							cal_flag[ib][id]=false;
							continue;
						}


						Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[ib][id][0], dr[ib][id][1], dr[ib][id][2], rly, grly);

						// 1E-7 is necessary in case of R is just on one grid
						const double position = distance[ib][id] / delta_r;
						this->iq[id] = static_cast<int>(position);
						this->x0[id] = position - static_cast<double>(iq[id]);
						this->x1[id] = 1.0 - x0[id];
						this->x2[id] = 2.0 - x0[id];
						this->x3[id] = 3.0 - x0[id];
						this->x12[id] = x1[id]*x2[id] / 6.0;
						this->x03[id] = x0[id]*x3[id] / 2.0;

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
									dphix[ib][id][iw] = dphiy[ib][id][iw] = dphiz[ib][id][iw] = 0.0;
								}
								else
								{
									pointer = &ORB.Phi[it].
										PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

									double Zty = pointer->zty;
									psir_ylm[ib][id][iw] = Zty * rly[idx_lm];
									dphix[ib][id][iw] = Zty * grly[idx_lm][0];
									dphiy[ib][id][iw] = Zty * grly[idx_lm][1];
									dphiz[ib][id][iw] = Zty * grly[idx_lm][2];
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

								dphix[ib][id][iw] = tmpdphi_rly * dr[ib][id][0]  + tmprl * grly[idx_lm][0];
								dphiy[ib][id][iw] = tmpdphi_rly * dr[ib][id][1]  + tmprl * grly[idx_lm][1];
								dphiz[ib][id][iw] = tmpdphi_rly * dr[ib][id][2]  + tmprl * grly[idx_lm][2];
							}
						}
					}// ib
				}//!id //finish loop of calc pre-info for each adjacent atom

                //get information about the vector potential A, and position r//
				int index=0;
				double rA[3]={0,0,0};     
				for(int ii=0; ii<pw.bx; ii++)
				{
					const double iii = (i*pw.bx + ii)/(double)pw.ncx;
					for(int jj=0; jj<pw.by; jj++)
					{
						const double jjj = (j*pw.by + jj)/(double)pw.ncy;
						for(int kk=0; kk<pw.bz; kk++)
						{
							const double kkk = (k*pw.bz + kk)/(double)pw.ncz;
							// A * r
							rA[0]=iii*latvec1[0]+jjj*latvec2[0]+kkk*latvec3[0];
							rA[1]=iii*latvec1[1]+jjj*latvec2[1]+kkk*latvec3[1];
							rA[2]=iii*latvec1[2]+jjj*latvec2[2]+kkk*latvec3[2];
                            Rbxyz[index][0]=rA[0]; 
                            Rbxyz[index][1]=rA[1]; 
                            Rbxyz[index][2]=rA[2]; 
							bfid.cal_A(Aldr3[index],rA);
							// A * A
							AAdr[index]=Aldr3[index][0]*Aldr3[index][0]+Aldr3[index][1]*Aldr3[index][1]+Aldr3[index][2]*Aldr3[index][2];   //add at 2012-01-06//                     
							index++;  //updated at 2012-01-05//
						}
					}
				}

                //calculate the matrix element using the pre-information above//
				for (int ia1=0; ia1<size; ia1++)
				{
					const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
					const int iat1 = gt.which_atom[mcell_index1];
					const int T1 = ucell.iat2it[iat1];
					const int I1 = ucell.iat2ia[iat1];
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
					Atom *atom1 = &ucell.atoms[T1];
					//get vector A of atom 1.
                    A_atom1[0] = bfid.A_of_Atom[iat1][0];
                    A_atom1[1] = bfid.A_of_Atom[iat1][1];
                    A_atom1[2] = bfid.A_of_Atom[iat1][2];
					for (int ia2=0; ia2<size; ia2++)
					{
						const int mcell_index2 = gt.bcell_start[grid_index] + ia2;
						const int iat2 = gt.which_atom[mcell_index2];
						const int T2 = ucell.iat2it[iat2];
						const int I2 = ucell.iat2ia[iat2];
						const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
						const int iw2_lo = gt.trace_lo[start2];
						Atom *atom2 = &ucell.atoms[T2];
					    //get vector A of atom 2.
						A_atom2[0] = bfid.A_of_Atom[iat2][0];
						A_atom2[1] = bfid.A_of_Atom[iat2][1];
						A_atom2[2] = bfid.A_of_Atom[iat2][2];
						//get A2-A1, zhiyuan add 2012-01-05//
						A2_ms_A1[0]=A_atom2[0]-A_atom1[0];
						A2_ms_A1[1]=A_atom2[1]-A_atom1[1];
						A2_ms_A1[2]=A_atom2[2]-A_atom1[2];
						for(int ib=0; ib<pw.bxyz; ib++)
						{
							if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
							{
								double* psi1 = &psir_ylm[ib][ia1][0];
								double* dphi1x = dphix[ib][ia1];
								double* dphi1y = dphiy[ib][ia1];
								double* dphi1z = dphiz[ib][ia1];
							    complex<double> phase=0;  //the phase for overlap//
								phase=bfid.cal_phase(A2_ms_A1,Rbxyz[ib],-1);

								int iw1_lo = gt.trace_lo[start1];
								
								//cycle for the orbitals of atom1//
								for (int iw=0; iw< atom1->nw; ++iw, ++psi1, 
									  ++dphi1x,++dphi1y,++dphi1z)
								{
									complex<double>* pS = &GridS[iw1_lo][iw2_lo];
									complex<double>* pT = &Grid_T[iw1_lo][iw2_lo];
									complex<double>* pAP = &GridAP[iw1_lo][iw2_lo];
									complex<double>* pAA = &GridAA[iw1_lo][iw2_lo];

									double* psi2= &psir_ylm[ib][ia2][0];
									double* dphi2x = dphix[ib][ia2];
									double* dphi2y = dphiy[ib][ia2];
									double* dphi2z = dphiz[ib][ia2];

									double* phi2_end = psi2 + atom2->nw;

								    //cycle for the orbitals of atom2//
									for (; psi2<phi2_end; ++pS,++pT,++pAP,++pAA,
											++psi2,++dphi2x, ++dphi2y, ++dphi2z)
									{
										*pS += Add_S(phase,*psi1,*psi2);
										*pT += Add_T(phase,A_atom1,A_atom2,*dphi1x,*dphi1y,*dphi1z,*dphi2x,*dphi2y,*dphi2z,*psi1,*psi2);
										*pAP += Add_AP(phase,Aldr3[ib],A_atom2,*dphi2x,*dphi2y,*dphi2z,*psi1,*psi2);
										*pAA += Add_AA(phase,AAdr[ib],*psi1,*psi2);
									}
									iw1_lo++;
								}//iw1
							}// cal_flag
						}//ib
					}//!ia
				}// ia1
			}// k
		}// j
	}// i
//===========================big cycle ends============================================

    //delete arrays//
	if(max_size!=0)
	{
		for(int i=0; i<pw.bxyz; i++)
		{
			for(int j=0; j<max_size; j++)
			{
				delete[] dr[i][j];
				delete[] psir_ylm[i][j];
				delete[] dphix[i][j];
				delete[] dphiy[i][j];
				delete[] dphiz[i][j];
			}
			delete[] dr[i];
			delete[] distance[i];
			delete[] psir_ylm[i];
			delete[] cal_flag[i];
			delete[] dphix[i];
			delete[] dphiy[i];
			delete[] dphiz[i];
		}
		delete[] dr;
		delete[] distance;
		delete[] psir_ylm;
		delete[] cal_flag;
		delete[] dphix;
		delete[] dphiy;
		delete[] dphiz;
	}

	//delete Aldr3, Rbxyz and AAdr//
	for(int i=0;i<pw.bxyz;i++)
	{
		delete [] Aldr3[i];
		delete [] Rbxyz[i];
	}
	delete [] Aldr3;
	delete [] Rbxyz;
	delete [] AAdr;

    //parallel to all processors, and upload to S,H(temporary) matrixes//
	for (int i=0; i<NLOCAL; i++)
	{
		complex<double> *tmpS = new complex<double>[NLOCAL];
		complex<double> *tmpT = new complex<double>[NLOCAL];
		complex<double> *tmpAP = new complex<double>[NLOCAL];
		complex<double> *tmpAA = new complex<double>[NLOCAL];

		const int mu = gt.trace_lo[i];
		if(mu >=0)
		{
			for (int j=0; j<NLOCAL; j++)
			{
				const int nu = gt.trace_lo[j];
				if(nu>=0)
				{
					tmpS[j] = GridS[mu][nu];
					tmpT[j] = Grid_T[mu][nu];
					tmpAP[j] = GridAP[mu][nu];
					tmpAA[j] = GridAA[mu][nu];
				}
			}
		}

        Parallel_Reduce::reduce_complex_double_pool( tmpS, NLOCAL); 
        Parallel_Reduce::reduce_complex_double_pool( tmpT, NLOCAL); 
        Parallel_Reduce::reduce_complex_double_pool( tmpAP, NLOCAL); 
        Parallel_Reduce::reduce_complex_double_pool( tmpAA, NLOCAL); 

		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			complex<double> SS=tmpS[j]; 
			complex<double> TT=tmpT[j]; 
			complex<double> AP=tmpAP[j]; 
			complex<double> AA=tmpAA[j]; 
			
			if(type=='S')
			{
				LM.set_HSk(i,j,SS,'S');
			}
			else if(type=='T')
			{
				LM.set_HSk(i,j,TT,'T');
			}
			else if(type=='A')
			{
				LM.set_HSk(i,j,AP,'T');
				LM.set_HSk(i,j,AA,'T');
			}
		}

		delete[] tmpS;
		delete[] tmpT;
		delete[] tmpAP;
		delete[] tmpAA;
	}
	
	//delete GridS, Grid_T, GridAP, GridAA//
	if(gt.lgd>0)
	{
		for(int i=0; i<gt.lgd; i++)
		{
			delete[] GridS[i];
			delete[] Grid_T[i];
			delete[] GridAP[i];
			delete[] GridAA[i];
		}
		delete[] GridS;
		delete[] Grid_T;
		delete[] GridAP;
		delete[] GridAA;
	}
	
	timer::tick("Gint_Gamma","gamma_S_T_AP");
	return;
}

//sub-function for calculating S//
complex<double> Gint_Gamma::Add_S(complex<double> phase, double psi1, double psi2)
{
	return (phase*psi1*psi2*this->vfactor);
}

//sub-function for calculating T//
complex<double> Gint_Gamma::Add_T(complex<double> phase, double A1[3], double A2[3], 
		double dphi1x,double dphi1y,double dphi1z,
		double dphi2x,double dphi2y,double dphi2z, 
		double psi1, double psi2)
{
	double part1 = 2.0/bfid.c/bfid.c*(A1[0]*A2[0]+A1[1]*A2[1]+A1[2]*A2[2])*psi1*psi2;
	double part2 = dphi1x*dphi2x+dphi1y*dphi2y+dphi1z*dphi2z; 
	double part3 = -1.414213562/bfid.c*(A2[0]*dphi1x+A2[1]*dphi1y+A2[2]*dphi1z)*psi2;
    double part4 =  1.414213562/bfid.c*(A1[0]*dphi2x+A1[1]*dphi2y+A1[2]*dphi2z)*psi1;
	complex<double> sum1= complex<double> ((part1+part2),(part3+part4));
	complex<double> sum = phase*sum1*this->vfactor;
	return(sum);
}

//sub-function for calculating A*P//
complex<double> Gint_Gamma::Add_AP(complex<double> phase, double Aldr3[3], double A2[3], 
                                   double dphi2x,double dphi2y,double dphi2z,
								   double psi1, double psi2)
{
    double part1 = (Aldr3[0]*dphi2x+Aldr3[1]*dphi2y+Aldr3[2]*dphi2z);
    double part2 = -1.414213562/bfid.c*(Aldr3[0]*A2[0]+Aldr3[1]*A2[1]+Aldr3[2]*A2[2])*psi2;
	complex<double> sum1 = complex<double> (part1,part2);
	complex<double> sum2 = -2.0*1.414213562/bfid.c*phase*sum1*psi1*this->vfactor;
	complex<double> sum = (complex<double> (0,1.0))*sum2;
	return(sum);
}

//sub-function for calculating A*A//
complex<double> Gint_Gamma::Add_AA(complex<double> phase, double AAdr,
                                   double psi1, double psi2)
{
	return (phase*AAdr*2.0/bfid.c/bfid.c*psi1*psi2*this->vfactor);
}
