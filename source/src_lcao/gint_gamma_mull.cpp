#include "gint_gamma.h"
#include "grid_technique.h"
//#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"

#include "global_fp.h" // mohan add 2021-01-30
#include "../module_base/ylm.h"

void Gint_Gamma::cal_mulliken(double** mulliken)
{
    TITLE("Gint_Gamma","cal_mulliken");
    timer::tick("Gint_Gamma","cal_mulliken");

    this->save_atoms_on_grid(GridT);
    this->gamma_mulliken(mulliken);

    timer::tick("Gint_Gamma","cal_mulliken");
    return;
}

// this subroutine lies in the heart of LCAO algorithms.
// so it should be done very efficiently, very carefully.
// I might repeat again to emphasize this: need to optimize
// this code very efficiently, very carefully.
void Gint_Gamma::gamma_mulliken(double** mulliken)
{
    TITLE("Grid_Integral","gamma_charge");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;

	// allocate 1
	int nnnmax=0;
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax = max(nnnmax, nnn[T]);
	}

	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;	
	bool** cal_flag;
	if(max_size!=0) 
	{
		dr = new double**[pw.bxyz];
		distance = new double*[pw.bxyz];
		psir_ylm = new double**[pw.bxyz];
		cal_flag = new bool*[pw.bxyz];

		for(int i=0; i<pw.bxyz; i++)
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

	double mt[3]={0,0,0};
	double *vldr3 = new double[pw.bxyz];
	double v1 = 0.0;
	int* vindex=new int[pw.bxyz];
	ZEROS(vldr3, pw.bxyz);
	ZEROS(vindex, pw.bxyz);
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

                // get the value: how many atoms has orbital value on this grid.
                const int size = GridT.how_many_atoms[ this->grid_index ];
				if(size==0) continue;

				// (1) initialized the phi * Ylm.
                for (int id=0; id<size; id++)
                {
                    // there are two parameters we want to know here:
                    // in which bigcell of the meshball the atom in?
                    // what's the cartesian coordinate of the bigcell?
                    const int mcell_index = GridT.bcell_start[grid_index] + id;
                    const int imcell = GridT.which_bigcell[mcell_index];

                    int iat = GridT.which_atom[mcell_index];

                    const int it = ucell.iat2it[ iat ];
                    const int ia = ucell.iat2ia[ iat ];

                    // meshball_positions should be the bigcell position in meshball
                    // to the center of meshball.
                    // calculated in cartesian coordinates
                    // the vector from the grid which is now being operated to the atom position.
                    // in meshball language, is the vector from imcell to the center cel, plus
                    // tau_in_bigcell.
					mt[0] = GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
					mt[1] = GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
					mt[2] = GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

					for(int ib=0; ib<pw.bxyz; ib++)
					{
						// meshcell_pos: z is the fastest
						dr[ib][id][0] = GridT.meshcell_pos[ib][0] + mt[0]; 
						dr[ib][id][1] = GridT.meshcell_pos[ib][1] + mt[1]; 
						dr[ib][id][2] = GridT.meshcell_pos[ib][2] + mt[2]; 	

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
						
						std::vector<double> ylma;
						//if(distance[id] > GridT.orbital_rmax) continue;
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
						/*
						   this->iq[id] = static_cast<int>(position);
						   this->x0[id] = position - static_cast<double>(iq[id]);
						   this->x1[id] = 1.0 - x0[id];
						   this->x2[id] = 2.0 - x0[id];
						   this->x3[id] = 3.0 - x0[id];
						   this->x12[id] = x1[id]*x2[id] / 6.0;
						   this->x03[id] = x0[id]*x3[id] / 2.0;
						 */
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

						//		int ip = this->iq[id];
						//		double A = ip+1.0-position/delta_r;
						//		double B = 1.0-A;
						//		double coef1 = (A*A*A-A)/6.0*delta_r*delta_r;
						//		double coef2 = (B*B*B-B)/6.0*delta_r*delta_r;

						Atom* atom1 = &ucell.atoms[it];
						for (int iw=0; iw< atom1->nw; iw++)
						{
							if ( atom1->iw2_new[iw] )
							{
								pointer = &ORB.Phi[it].PhiLN(
										atom1->iw2l[iw],
										atom1->iw2n[iw]);
								phi = c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
									+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
							}
							psir_ylm[ib][id][iw] = phi * ylma[atom1->iw2_ylm[iw]];
							//psir_ylm[ib][id][iw] = 1;//for test
						}
					}// end ib
                }// end id

                for (int ia1=0; ia1<size; ia1++)
                {
                    const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
					const int T1 = ucell.iat2it[ GridT.which_atom[mcell_index1] ];
					Atom *atom1 = &ucell.atoms[T1];
					const int I1 = ucell.iat2ia[ GridT.which_atom[mcell_index1] ];
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
						const int mcell_index2 = GridT.bcell_start[grid_index] + ia2;
						const int T2 = ucell.iat2it[ GridT.which_atom[mcell_index2]];

						// only do half part of matrix(including diago part)
						// for T2 > T1, we done all elements, for T2 == T1,
						// we done about half.
						if (T2 >= T1)
						{
							Atom *atom2 = &ucell.atoms[T2];
							const int I2 = ucell.iat2ia[ GridT.which_atom[mcell_index2]];
							const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

							for(int is=0; is<NSPIN; is++)
							{
								double *rhop = CHR.rho[is];
								for (int ib=0; ib<pw.bxyz; ib++)
								{
									if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
									{
										int iw1_lo = GridT.trace_lo[start1];
										int iw1_all = start1;
										double* psi1 = psir_ylm[ib][ia1];
										double* psi2 = psir_ylm[ib][ia2];

										// how many orbitals in this type: SZ or DZP or TZP...
										for (int iw=0; iw< atom1->nw; iw++, ++iw1_lo, ++iw1_all)
										{
											v1=psi1[iw]+psi1[iw];
											int iw2_lo = GridT.trace_lo[start2];
											double *DMp = &LOC.DM[is][iw1_lo][iw2_lo];
											double *psi2p = psi2;
											double *psi2p_end = psi2 + atom2->nw;

											double tmp = 0.0;
											for (; psi2p < psi2p_end; ++iw2_lo, ++psi2p, ++DMp)
											{
												if ( iw1_lo > iw2_lo)
												{
													continue;
												}
												// for diago part, the charge density be accumulated once.
												// for off-diago part, the charge density be accumulated twice.
												// (easy to understand, right? because we only calculate half
												// of the matrix).

												double tmp1 = v1 * psi2p[0];
												double tmp2 = tmp1 * DMp[0];
												if (iw1_lo<iw2_lo)
												{
													tmp += tmp2;
												}
												else
												{
													tmp += tmp2/2.0;
												}
											}//iw2
											mulliken[is][iw1_all] += tmp;
										}//iw
									}// cal_flag
								}//ib
							}
                        }//T
                    }// ia2
                }// ia1

            }// k
        }// j
    }// i

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
			delete[] cal_flag[i];
		}
		delete[] dr;
		delete[] distance;
		delete[] psir_ylm;
		delete[] cal_flag;
	}

	for(int is=0; is<NSPIN; ++is)
	{
		Parallel_Reduce::reduce_double_all( mulliken[is], NLOCAL );
	}

    return;
}

