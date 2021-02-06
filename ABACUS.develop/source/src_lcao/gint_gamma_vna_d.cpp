#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"

#include "global_fp.h" // mohan add 2021-01-30

void Gint_Gamma::cal_vna_d(const Grid_Technique &gt, const double* vrs1,
const char &matrix_type)
{
    TITLE("Gint_Gamma","cal_vna_d");
    timer::tick("Gint_Gamma","cal_vna_d");

    this->save_atoms_on_grid(gt);
	this->gamma_vna_d(gt, matrix_type);

    timer::tick("Gint_Gamma","cal_vna_d");
    return;
}

void Gint_Gamma::gamma_vna_d(const Grid_Technique &gt, const char &matrix_type)
{
    TITLE("Grid_Integral","gamma_vna_d");
//	cout << " gamma_vna_d ()" << endl;

	double** GridVlocal;
	bool perform_gint = true;

	const int lgd_now = gt.lgd;
	if(lgd_now > 0)
	{
		GridVlocal = new double*[lgd_now];
		for (int i=0; i<lgd_now; i++)
		{
			GridVlocal[i] = new double[lgd_now];
			ZEROS(GridVlocal[i], lgd_now);
		}
		Memory::record("Gint_Gamma","GridVlocal",lgd_now*lgd_now,"double");
	}
	else if(lgd_now <= 0)
	{
		perform_gint = false;
	}

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;
	if( delta_r != ORB.Vna[0].dr )
	{
		cout << " delta_r = " << delta_r << endl;
		cout << " dr = " << ORB.Vna[0].dr << endl;
		WARNING_QUIT("Gint_Gamma::gamma_vna_d","delta_r != ORB.Vna[0].dr");
	}

	// allocate 1
	int nnnmax=0;
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax = max(nnnmax, nnn[T]);
	}
	
	const int bxyz = gt.bxyz;
	const int bx = gt.bx;
	const int by = gt.by;
	const int bz = gt.bz;

	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;	
	bool** cal_flag;
	if(max_size!=0) 
	{
		dr = new double**[bxyz];
		distance = new double*[bxyz];
		psir_ylm = new double**[bxyz];
		cal_flag = new bool*[bxyz];

		for(int i=0; i<bxyz; i++)
		{
			dr[i] = new double*[max_size];
			distance[i] = new double[max_size];
			psir_ylm[i] = new double*[max_size];
			cal_flag[i] = new bool[max_size];

			ZEROS(distance[i], max_size);
			ZEROS(cal_flag[i], max_size);

			for(int j=0; j<max_size; j++) 
			{
				dr[i][j] = new double[3];
				psir_ylm[i][j] = new double[ucell.nwmax];
				ZEROS(dr[i][j],3);
				ZEROS(psir_ylm[i][j],ucell.nwmax);
			}
		}
	}

	double* ylma = new double[nnnmax]; // Ylm for each atom: [bxyz, nnnmax]
	ZEROS(ylma, nnnmax);

	double mt[3]={0,0,0};
	double *vldr3 = new double[bxyz];
	double v1=0.0;
	ZEROS(vldr3, bxyz);
	double phi=0.0;

	
	double* vna3d = new double[bxyz];
	ZEROS(vna3d, bxyz);

	const int nbx = gt.nbx;
	const int nby = gt.nby;
	const int nbz_start = gt.nbzp_start;
	const int nbz = gt.nbzp;
	const int ncy = gt.ncy;

	//ofstream ofs("vna_d.txt");

    for (int i=0; i<nbx; ++i)
    {
		if( !perform_gint ) continue;
        for (int j=0; j<nby; ++j)
        {
            for (int k=nbz_start; k<nbz_start+nbz; ++k) // FFT grid
            {
                this->grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;

                // get the value: how many atoms has orbital value on this grid.
                const int size = gt.how_many_atoms[ this->grid_index ];
				if(size==0) continue;

				ZEROS(vna3d, bxyz);

				// (1) initialized the phi * Ylm.
                for (int id=0; id<size; ++id)
                {
                    // there are two parameters we want to know here:
                    // in which bigcell of the meshball the atom in?
                    // what's the cartesian coordinate of the bigcell?
                    const int mcell_index = gt.bcell_start[grid_index] + id;
                    const int imcell = gt.which_bigcell[mcell_index];

                    int iat = gt.which_atom[mcell_index];

                    const int it = ucell.iat2it[ iat ];
                    const int ia = ucell.iat2ia[ iat ];

                    // meshball_positions should be the bigcell position in meshball
                    // to the center of meshball.
                    // calculated in cartesian coordinates
                    // the vector from the grid which is now being operated to the atom position.
                    // in meshball language, is the vector from imcell to the center cel, plus
                    // tau_in_bigcell.
					mt[0] = gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0];
					mt[1] = gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1];
					mt[2] = gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2];

					for(int ib=0; ib<bxyz; ib++)
					{
						// meshcell_pos: z is the fastest
						dr[ib][id][0] = gt.meshcell_pos[ib][0] + mt[0]; 
						dr[ib][id][1] = gt.meshcell_pos[ib][1] + mt[1]; 
						dr[ib][id][2] = gt.meshcell_pos[ib][2] + mt[2]; 	

						distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] + 
						dr[ib][id][1]*dr[ib][id][1] 
						+ dr[ib][id][2]*dr[ib][id][2]);

						if(distance[ib][id] > ORB.Vna[it].rcut &&
							distance[ib][id] > ORB.Phi[it].getRcut())
						{
							cal_flag[ib][id]=false;
							continue;
						}

						//if(distance[id] > gt.orbital_rmax) continue;
						//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
						if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
						
						Ylm::sph_harm (	ucell.atoms[it].nwl,
								dr[ib][id][0] / distance[ib][id],
								dr[ib][id][1] / distance[ib][id],
								dr[ib][id][2] / distance[ib][id],
								ylma);
						// these parameters are about interpolation
						// because once we know the distance from atom to grid point,
						// we can get the parameters we need to do interpolation and
						// store them first!! these can save a lot of effort.
						const double position = distance[ib][id] / delta_r;
						this->iq[id] = static_cast<int>(position);
						this->x0[id] = position - static_cast<double>(iq[id]);
						this->x1[id] = 1.0 - x0[id];
						this->x2[id] = 2.0 - x0[id];
						this->x3[id] = 3.0 - x0[id];
						this->x12[id] = x1[id]*x2[id] / 6.0;
						this->x03[id] = x0[id]*x3[id] / 2.0;

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
						
						double ccc;
						if(distance[ib][id] <= ORB.Vna[it].rcut)
						{
							ccc = (x12[id]*(ORB.Vna[it].vna_u[iq[id]]*x3[id]
							+ORB.Vna[it].vna_u[iq[id]+3]*x0[id])
							+ x03[id]*(ORB.Vna[it].vna_u[iq[id]+1]*x2[id]
							-ORB.Vna[it].vna_u[iq[id]+2]*x1[id]));
						//	cout << " ccc=" << ccc << endl;
						//	int ok; cin >> ok;
							vna3d[ib] += ccc;
						}

//						if(i==0 && j==7 &&k==13 && ib==6)
//						{
//							ofs << " dis=" << distance[ib][id] << " ccc=" << ccc << endl;
//							ofs << " 1=" << ORB.Vna[it].vna_u[iq[id]] << endl;
//							ofs << " 2=" << ORB.Vna[it].vna_u[iq[id]+1] << endl;
//							ofs << " 3=" << ORB.Vna[it].vna_u[iq[id]+2] << endl;
//							ofs << " 4=" << ORB.Vna[it].vna_u[iq[id]+3] << endl;
//						}

						//		int ip = this->iq[id];
						//		double A = ip+1.0-position/delta_r;
						//		double B = 1.0-A;
						//		double coef1 = (A*A*A-A)/6.0*delta_r*delta_r;
						//		double coef2 = (B*B*B-B)/6.0*delta_r*delta_r;


						if(distance[ib][id] <= ORB.Phi[it].getRcut()) 
						{
							cal_flag[ib][id]=true;
						}
						else 
						{
							cal_flag[ib][id]=false;
							continue;
						}
	

						Atom* atom1 = &ucell.atoms[it];
						for (int iw=0; iw< atom1->nw; iw++)
						{
							if ( atom1->iw2_new[iw] )
							{
								pointer = &ORB.Phi[it].PhiLN(
										atom1->iw2l[iw],
										atom1->iw2n[iw]);
								phi = c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
									+ c3*pointer->psi_uniform[ip+1]+c4*pointer->dpsi_uniform[ip+1];
							}
							psir_ylm[ib][id][iw] = phi * ylma[atom1->iw2_ylm[iw]];
							//psir_ylm[ib][id][iw] = 1;//for test
						}
					}// end ib
                }// end id


				// extract the local potentials.
				for(int ib=0; ib<bxyz; ib++)
				{
				//	vna3d[ib] = ib;//tmp by mohan
					vldr3[ib] = vna3d[ib]  * this->vfactor;

//					if(i==0 && j==7 &&k==13 && ib==6)
//					{
//						ofs << " i=" << i << " j=" << j << " k=" << k << " ib=" << ib << " vna3d=" << vna3d[ib] << endl;
//						for(int id=0; id<size; id++)
//						{
//							ofs << " " << id << " " << distance[ib][id] << endl;
//						}
//					}

					// vldr3[ib] = 1.0e-5; // for test
					// vldr3[bindex] = this->vfactor; // for checking overlap S
				}

                for (int ia1=0; ia1<size; ia1++)
                {
                    const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
					const int iat1 = gt.which_atom[mcell_index1];
					const int T1 = ucell.iat2it[iat1];
					Atom *atom1 = &ucell.atoms[T1];
					const int I1 = ucell.iat2ia[iat1];
					// get the start index of local orbitals.
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

					// call to get real spherical harmonic values according to
					// a particular number: (lmax+1)^2 and vectors between
					// atom and grid point(we need the direction), the output
					// are put in the array: ylm1.
					//Ylm::get_ylm_real(this->nnn[T1], this->dr[ia1], this->ylm1);

					// attention! assume all rcut are same for this atom type now.
					//if (distance[ia1] > ORB.Phi[T1].getRcut())continue;

					//for(int ia2=ia1; ia2<size; ia2++)
					for (int ia2=0; ia2<size; ia2++)
					{
						const int mcell_index2 = gt.bcell_start[grid_index] + ia2;
						const int T2 = ucell.iat2it[ gt.which_atom[mcell_index2]];

						// only do half part of matrix(including diago part)
						// for T2 > T1, we done all elements, for T2 == T1,
						// we done about half.
						if (T2 >= T1)
						{
							Atom *atom2 = &ucell.atoms[T2];
							const int I2 = ucell.iat2ia[ gt.which_atom[mcell_index2]];
							const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

							for (int ib=0; ib<bxyz; ib++)
							{
								if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
								{
									double* psi1 = psir_ylm[ib][ia1];
									double* psi2 = psir_ylm[ib][ia2];
									int iw1_lo = gt.trace_lo[start1];

									// how many orbitals in this type: SZ or DZP or TZP...
									for (int iw=0; iw< atom1->nw; iw++, ++iw1_lo)
									{
										v1=psi1[iw] * vldr3[ib];
										int iw2_lo = gt.trace_lo[start2];
										double* result = &GridVlocal[iw1_lo][iw2_lo];
										double* psi2p = psi2;
										double* psi2p_end = psi2p + atom2->nw;
										for (;psi2p<psi2p_end; ++psi2p, ++result, ++iw2_lo)
										{
											// we only need to calculate half of the matrix
											// (including diago part). what we need to calculate
											// are the terms satisfy the condition:
											// iw1_all <= iw2_all.
											// same for iw1_lo <= iw2_lo;
											// the matrix is always symmetry.
											if ( iw1_lo > iw2_lo)
											{
												continue;
											}

											result[0] += v1 * psi2p[0];

											// note 1
											// why I increase both iw2_lo and iw2_all:
											// Because iw2_lo is used to save element.
											// And iw2_all is used to judge if we need
											// to calculate the element according to
											// "iw1_all > iw2_all" condition.
											// note 2
											// in fact we don't need to do this.
											// because GridVlocal is a symmetry matrix.
											// whatever the division is , the order
											// is same between iw2_lo,iw1_lo and
											// iw2_all, iw1_all
										}//iw2
                            		}//iw
								}// cal_flag
							}//ib
                        }//T
                    }// ia2
                }// ia1

            }// k
        }// j
    }// i

	delete[] ylma;
	delete[] vldr3;
	delete[] vna3d;
	
	if(max_size!=0) 
	{
		for(int i=0; i<bxyz; i++)
		{
			for(int j=0; j<max_size; j++) 
			{
				delete[] dr[i][j];
				delete[] psir_ylm[i][j];
			}
			delete[] dr[i];
			delete[] distance[i];
			delete[] psir_ylm[i];
			delete[] cal_flag[i];
		}
		delete[] dr;
		delete[] distance;
		delete[] psir_ylm;
		delete[] cal_flag;
	}

    {
        double* tmp;
        for (int i=0; i<NLOCAL; i++)
        {
            tmp = new double[NLOCAL];
            ZEROS(tmp, NLOCAL);
            const int mu = gt.trace_lo[i];
            if (mu >= 0)
            {
                for (int j=0; j<NLOCAL; j++)
                {
                    const int nu = gt.trace_lo[j];
                    if (nu >=0)
                    {
                        if (mu <= nu)
                        {
                            tmp[j] = GridVlocal[mu][nu];
                        }
                        else
                        {
							//-------------------------------
							// origin:
                            // tmp[i] = GridVlocal[nu][mu];
                            // mohan fix bug 
							// 2011-01-13
							//-------------------------------
							tmp[j] = GridVlocal[nu][mu];
                        }
                    }
                }
            }
            Parallel_Reduce::reduce_double_pool( tmp, NLOCAL );
            for (int j=0; j<NLOCAL; j++)
            {
                if (!ParaO.in_this_processor(i,j))
                {
                    continue;
                }
			
				// mohan update 2011-04-15	
				//change L to T, which is two different matrix.
				// L : local
				// T : T, Vnl, Vna
               	LM.set_HSgamma(i,j,tmp[j],matrix_type);
            }
            delete[] tmp;
        }
    }

	// mohan update 2010-09-07
	if(gt.lgd>0)
	{
		for(int i=0; i<gt.lgd; i++)
		{
			delete[] GridVlocal[i];
		}
		delete[] GridVlocal;
	}


    return;
}

