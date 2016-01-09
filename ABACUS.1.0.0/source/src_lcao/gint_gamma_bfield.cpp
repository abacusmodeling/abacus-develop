#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "bfield.h"

void Gint_Gamma::cal_vnl_B(complex<double> **Tab)
{
	TITLE("Gint_Gamma","cal_vnl_B");
	timer::tick("Gint_Gamma","cal_vnl_B");
	
	this->job = cal_vnlb;
	this->save_atoms_on_grid(GridT);
	this->grid_integration_vnl(Tab);

	timer::tick("Gint_Gamma","cal_vnl_B");
	return;
}

void Gint_Gamma::grid_integration_vnl(complex<double> **Tab)
{
    TITLE("Grid_Integral","grid_integration_vnl");
	cout << " grid integration for vnl() " << endl;
	cout << " vfactor=" << vfactor << endl;
	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	
	const int lgd_now= GridT.lgd;
	const int lgd_beta=GridT.lgbeta;

	complex<double> **table;
	if(lgd_now > 0 && lgd_beta > 0)
	{
		table = new complex<double>*[lgd_now];
		for (int i=0; i<lgd_now; i++)
		{
			table[i] = new complex<double> [lgd_beta];
			ZEROS(table[i], lgd_beta);
		}
	}
    
	const double delta_r = ORB.dr_uniform;
	const Numerical_Orbital_Lm* pa;
	const Numerical_Nonlocal_Lm* pb; 

	int nnnmax=0;
	int nhmax=0;//mohan add 2012-04-13
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax = max(nnnmax, nnn[T]);
		nhmax = max(nhmax, ucell.atoms[T].nh);
	}
	cout << " nhmax = " << nhmax << endl;

	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;
	double*** betar_ylm;
	bool** cal_flag;

	if(max_size!=0)
	{
		dr = new double**[pw.bxyz];
		distance = new double*[pw.bxyz];
		psir_ylm = new double**[pw.bxyz];
		betar_ylm = new double**[pw.bxyz];//mohan add 2012-04-13
		cal_flag = new bool*[pw.bxyz];
		for(int i=0; i<pw.bxyz; i++)
		{
			dr[i] = new double*[max_size];
			distance[i] = new double[max_size];
			psir_ylm[i] = new double*[max_size];
			betar_ylm[i] = new double*[max_size];//mohan add 2012-04-13
			cal_flag[i] = new bool[max_size];
			ZEROS(distance[i], max_size);
			ZEROS(cal_flag[i], max_size);
			for(int j=0; j<max_size; j++)
			{
				dr[i][j] = new double[3];
				psir_ylm[i][j] = new double[ucell.nwmax];
				betar_ylm[i][j] = new double[nhmax];//mohan add 2012-04-13
				ZEROS(dr[i][j],3);
				ZEROS(psir_ylm[i][j],ucell.nwmax);
				ZEROS(betar_ylm[i][j],nhmax);
			}
		}
	}

	double* ylma = new double[nnnmax];
	ZEROS(ylma, nnnmax);

	double mt[3]={0,0,0};
	double *vldr3 = new double[pw.bxyz];
	double v1=0.0;
	int* vindex=new int[pw.bxyz];
	ZEROS(vldr3, pw.bxyz);
	ZEROS(vindex, pw.bxyz);


	// PHASE ON GRID
    // This is to store the position of grids, zhiyuan add 2012-02-11//
    double **Rbxyz = new double*[pw.bxyz]; // bxyz*3//
    for(int i=0;i<pw.bxyz;i++)
    {
        Rbxyz[i]=new double[3];
        ZEROS(Rbxyz[i], 3);
    }	
    //sun zhiyuan add definition 2011-08-14
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
	

	double phi=0.0;

	const int nbx = GridT.nbx;
	const int nby = GridT.nby;
	const int nbz_start = GridT.nbzp_start;
	const int nbz = GridT.nbzp;

    for (int i=0; i<nbx; i++)
    {
        for (int j=0; j<nby; j++)
        {
            for (int k=nbz_start; k<nbz_start+nbz; k++) // FFT grid
            {
                this->grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
                const int size = GridT.how_many_atoms[ this->grid_index ];
				if(size==0) continue;

                for (int id=0; id<size; id++)
                {
                    const int mcell_index = GridT.bcell_start[grid_index] + id;
                    const int imcell = GridT.which_bigcell[mcell_index];
                    int iat = GridT.which_atom[mcell_index];
                    const int it = ucell.iat2it[ iat ];
                    const int ia = ucell.iat2ia[ iat ];
                    mt[0] = GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
                    mt[1] = GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
                    mt[2] = GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

					for(int ib=0; ib<pw.bxyz; ++ib)
					{
						dr[ib][id][0] = GridT.meshcell_pos[ib][0] + mt[0];
						dr[ib][id][1] = GridT.meshcell_pos[ib][1] + mt[1];
						dr[ib][id][2] = GridT.meshcell_pos[ib][2] + mt[2];

						distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] + dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);
						if(distance[ib][id] <= ORB.Phi[it].getRcut())
						{
							cal_flag[ib][id]=true;
						}
						else
						{
							cal_flag[ib][id]=false;
							continue;
						}
						if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
						Ylm::sph_harm ( ucell.atoms[it].nwl,
						dr[ib][id][0] / distance[ib][id],
						dr[ib][id][1] / distance[ib][id],
						dr[ib][id][2] / distance[ib][id],
						ylma);
						const double position = distance[ib][id] / delta_r;
						int ip;
						double dx, dx2, dx3;
						double c1, c2, c3, c4;
						ip = static_cast<int>(position);
						dx = position - ip;
						dx2 = dx * dx;
						dx3 = dx2 * dx;
						c3 = 3.0*dx2-2.0*dx3;
						c1 = 1.0-c3;
						c2 = (dx-2.0*dx2+dx3)*delta_r;
						c4 = (dx3-dx2)*delta_r;
						Atom* atom1 = &ucell.atoms[it];
						double tmp;
						for (int iw=0; iw< atom1->nw; iw++)
						{
							if ( atom1->iw2_new[iw] )
							{
								pa = &ORB.Phi[it].PhiLN(
										atom1->iw2l[iw],
										atom1->iw2n[iw]);
								tmp = c1*pa->psi_uniform[ip]+c2*pa->dpsi_uniform[ip]
									+ c3*pa->psi_uniform[ip+1] + c4*pa->dpsi_uniform[ip+1];
							}
							psir_ylm[ib][id][iw] = tmp * ylma[atom1->iw2_ylm[iw]];
						}

						// calculation of betar_ylm, new codes;
						int ibm = 0;
						for(int iproj=0; iproj< ORB.nproj[it]; iproj++)
						{
							pb = &ORB.Beta[it].Proj[iproj];
							if (distance[ib][id] >= pb->getRcut() )
							{
								tmp = 0.0;
							}
							else
							{
								if( ip+1 >= pb->nr_uniform )
								{
									cout << " iproj=" << iproj << endl;
									cout << " ip+1=" << ip+1 << endl;
									cout << " nr_uniform=" << pb->nr_uniform << endl;
									WARNING_QUIT("grid_integration_vnl","ip+1 >= pb->nr_uniform");			
								} 

								tmp = c1*pb->beta_uniform[ip]+c2*pb->dbeta_uniform[ip]
									+ c3*pb->beta_uniform[ip+1] + c4*pb->dbeta_uniform[ip+1];
							}

							const int L = pb->getL();
							for(int m=0; m<2*L+1; m++)
							{
								assert( nnn[it] > ORB.ib2_ylm(it,ibm) );
								betar_ylm[ib][id][ibm] = tmp * ylma[ ORB.ib2_ylm(it,ibm) ]; 
								++ibm;
							}
						}// end iproj
					}// end ib
                }// end id


				// PHASE
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
                            rA[0]=iii*latvec1[0]+jjj*latvec2[0]+kkk*latvec3[0];
                            rA[1]=iii*latvec1[1]+jjj*latvec2[1]+kkk*latvec3[1];
                            rA[2]=iii*latvec1[2]+jjj*latvec2[2]+kkk*latvec3[2];
                            Rbxyz[index][0]=rA[0];
                            Rbxyz[index][1]=rA[1];
                            Rbxyz[index][2]=rA[2];
                            index++;  //updated at 2012-01-05//
                        }
                    }
                }

			

				for (int ia1=0; ia1<size; ia1++)
                {
                    const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
                    const int T1 = ucell.iat2it[ GridT.which_atom[mcell_index1] ];
                    Atom *atom1 = &ucell.atoms[T1];
                    const int I1 = ucell.iat2ia[ GridT.which_atom[mcell_index1] ];
                    const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
					const int iat1 = ucell.itia2iat(T1,I1);
                    
					for (int ia2=0; ia2<size; ia2++)
					{
						const int mcell_index2 = GridT.bcell_start[grid_index] + ia2;
						const int T2 = ucell.iat2it[ GridT.which_atom[mcell_index2]];
						
						Atom *atom2 = &ucell.atoms[T2];
						const int I2 = ucell.iat2ia[ GridT.which_atom[mcell_index2]];
						const int iat2 = ucell.itia2iat(T2,I2); 

						int start2 =ORB.itiaib2ib_all(T2, I2, 0);

						//get A2-A1, zhiyuan add 2012-01-05//
						double A2_ms_A1[3]={0,0,0}; // PHASE
						A2_ms_A1[0]= bfid.A_of_Atom[iat2][0]-bfid.A_of_Atom[iat1][0];
						A2_ms_A1[1]= bfid.A_of_Atom[iat2][1]-bfid.A_of_Atom[iat1][1];
						A2_ms_A1[2]= bfid.A_of_Atom[iat2][2]-bfid.A_of_Atom[iat1][2];

						for (int ib=0; ib<pw.bxyz; ++ib)
						{
							int iw_all = start1;
							if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
							{
								double* psi1 = psir_ylm[ib][ia1];
								double* beta2 = betar_ylm[ib][ia2];
								int iw1_lo = GridT.trace_lo[start1];

								complex<double> phase=ZERO;// PHASE
								phase=bfid.cal_phase(A2_ms_A1,Rbxyz[ib],-1);// PHASE
								for (int iw=0; iw< atom1->nw; iw++)
								{
									complex<double> v1 = phase*psi1[iw]*vfactor;//PHASE

									//double phase = Bfield::a_dot_b( this->dr[ia2], bfid.RxB[iat2] );
									//double e1 = std::cos(phase);
									//double e2 = std::sin(phase);
									//complex<double> v1e = v1 * complex<double>(e1, e2);
									int iproj_all = start2;
									const int iw1_lo = GridT.trace_lo[iw_all];
									for (int iproj=0; iproj< atom2->nh; iproj++)
									{
										assert(iw_all<NLOCAL);
										assert(iproj_all<ORB.nkb);
										//Tab[iw_all][ib_all] = Tab[iw_all][ib_all]+v1e*betar_ylm[ia2][ib];

										const int iw2_lo = GridT.trace_beta[iproj_all];
										assert(iw2_lo>=0);
										table[iw1_lo][iw2_lo] += v1 * beta2[iproj];

										//cout << " Tab=" << Tab[iw_all][ib_all] << " betar_ylm=" << betar_ylm[ia2][ib] << endl;
										//	int ok; cin >> ok;
										++iproj_all;
									}//iw2
									++iw_all;
								}//iw
							}//end cal_flag
						}// end ib
					}// ia2
                }// ia1
            }// k
        }// j
    }// i

    delete[] vindex;
    delete[] ylma;
    delete[] vldr3;

    if(max_size!=0)
    {
        for(int i=0; i<pw.bxyz; i++)
        {
            for(int j=0; j<max_size; j++)
            {
                delete[] dr[i][j];
                delete[] psir_ylm[i][j];
            }
            delete[] dr[i];
            delete[] distance[i];
            delete[] psir_ylm[i];
			delete[] betar_ylm[i];
            delete[] cal_flag[i];
        }
        delete[] dr;
        delete[] distance;
        delete[] psir_ylm;
		delete[] betar_ylm;
        delete[] cal_flag;
    }
	
	// reduce the table.
	int nkb=ORB.nkb;
	int index=0;
	complex<double> *tmp=new complex<double> [nkb];
	for(int i=0;i<NLOCAL;i++)
	{
		ZEROS(tmp,nkb);
		const int mu=GridT.trace_lo[i];
		if(mu>=0)
		{
			for(int j=0;j<nkb;j++)
			{
				const int nu=GridT.trace_beta[j];
				if(nu>=0)
				{
					tmp[j]=table[mu][nu];
				}
			}
		}
		Parallel_Reduce::reduce_complex_double_pool( tmp, nkb); //need to be updated//
		if(ParaO.trace_loc_row[i]>=0|| ParaO.trace_loc_col[i]>=0)
		{
			for(int j=0;j<nkb;j++)
			{
				Tab[index][j]=tmp[j];
			}
			++index;
		}
	}
	delete[] tmp;



	return;
}
