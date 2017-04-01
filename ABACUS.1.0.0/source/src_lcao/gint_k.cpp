#include "../src_pw/tools.h"
#include "gint_k.h"
#include "lcao_nnr.h"
#include "lcao_orbitals.h"
#include "grid_technique.h"
#include "ylm.h"
#include "../src_pw/global.h"

Gint_k::Gint_k()
{
	ik_now = 0;	
	pvpR_alloc_flag = false;
	pvnapR_alloc_flag = false;
	spin_now = -1; // for a start value, must not equal 1,2 or 4.
	reduced = true;// use reduced memory for H storage.
}

Gint_k::~Gint_k()
{
	
}

void Gint_k::reset_spin(const int &spin_now_in)
{
	this->spin_now = spin_now_in;
	return;
}



// this function is different from allocate_pvpR,
// this function is used once in an ionic iteration.
// the other function is used once in each electron iteration.
// nnrg may change when the ions move.
// be called in local_orbital_ions.cpp
void Gint_k::allocate_pvnapR(void)
{
	TITLE("Gint_k","allocate_pvnapR");

	assert( VNA>0 );

	if(this->pvnapR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::allocate_pvnapR","pvnapR has been allocated!");
	}

	assert(LNNR.nnrg!=0); //mohan update 2012-07-01
	this->pvnapR_reduced = new double[LNNR.nnrg];
	ZEROS( pvnapR_reduced, LNNR.nnrg);

	double mem = Memory::record("allocate_pvpR", "pvnapR_reduced", LNNR.nnrg , "double");
	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") ofs_running << " Memory of pvnapR : " << mem << " MB" << endl;
	if( mem > 800 )
	{
		ofs_warning << " memory for pvnapR = " << mem << endl;
		ofs_warning << " which is larger than 800 MB ! " << endl;
		WARNING_QUIT("Gint_k","allocate_pvpR");
	}

	this->pvnapR_alloc_flag = true;

	return;
}



void Gint_k::allocate_pvpR(void)
{
	TITLE("Gint_k","allocate_pvpR");

	if(this->pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::allocate_pvpR","pvpR has been allocated!");
	}

//	reduced = NURSE; 
	//xiaohui modify 2015-05-30
	//cout << " reduced algorithm for grid integration = " << reduced << endl;

	if(this->reduced)
	{
		// the number of matrix element <phi_0 | V | phi_R> is LNNR.nnrg.
		this->pvpR_reduced = new double[LNNR.nnrg];	
		ZEROS( pvpR_reduced, LNNR.nnrg);

		double mem = Memory::record("allocate_pvpR", "pvpR_reduced", LNNR.nnrg , "double");
	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
        if(OUT_LEVEL != "m") ofs_running << " Memory of pvpR : " << mem << " MB" << endl;
        if( mem > 800 )
        {
            ofs_warning << " memory for pvpR = " << mem << endl;
            ofs_warning << " which is larger than 800 MB ! " << endl;
			WARNING_QUIT("Gint_k","allocate_pvpR");
		}

	}
	else
	{
		double mem = Memory::record("allocate_pvpR", "pvpR", GridT.lgd * GridT.nutot
				* GridT.lgd * GridT.nutot , "double");
		//xiaohui add 'OUT_LEVEL' line, 2015-09-16
		if(OUT_LEVEL != "m") ofs_running << " Memory of pvpR : " << mem << " MB" << endl;
		if( mem > 800 )
		{
			ofs_warning << " memory for pvpR = " << mem << endl;
			ofs_warning << " which is larger than 800 MB ! " << endl;
			WARNING_QUIT("Gint_k","allocate_pvpR");
		}

		// output information
		//cout << " MEMORY OF pvpR       : " << Memory::record("allocate_pvpR", "pvpR", GridT.lgd * GridT.nutot 
		//* GridT.lgd * GridT.nutot , "cdouble") << " MB" << endl;

		//----------------------------------------------
		// allocate the complex matrix !!
		// nutot : total number of unitcells involved.
		// this may be very large, at least 
		// 3*3*3 = 27.
		//----------------------------------------------
		const int LDIM=GridT.lgd*GridT.nutot;
		this->pvpR_pool = new double[LDIM*LDIM];
		ZEROS(pvpR_pool, LDIM*LDIM);
		this->pvpR = new double*[LDIM];
		for(int i=0; i<LDIM; i++)
		{
			pvpR[i] = &pvpR_pool[i*LDIM];
		}
	}

	this->pvpR_alloc_flag = true;
	return;
}


void Gint_k::destroy_pvnapR(void)
{
	TITLE("Gint_k","destroy_pvnapR");
	
	if(!pvnapR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvnapR","<phi_0i | Vna | phi_Rj> matrix has not been allocated yet!");
	}
	
	delete[] pvnapR_reduced;

	this->pvnapR_alloc_flag = false;
	return;
}




void Gint_k::destroy_pvpR(void)
{
	TITLE("Gint_k","destroy_pvpR");
	
	if(!pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvpR","<phi_0i | V | phi_Rj> matrix has not been allocated yet!");
	}
	
	if(this->reduced)
	{
		delete[] pvpR_reduced;
	}	
	else
	{
		delete[] pvpR;
		delete[] pvpR_pool;
	}

	this->pvpR_alloc_flag = false;
	return;
}





// fold the <phi | vl |dphi(R)> * DM(R) to 
// calculate the force.
void Gint_k::folding_force(double** fvl_dphi,
	double* pvdpx, double* pvdpy, double* pvdpz)
{
	TITLE("Gint_k","folding_force");
	timer::tick("Gint_k","folding_force");

	//xiaohui modify 2013-12-17, test
//	assert(GridT.lgd > 0); //mohan add 2012-06-10

	// mohan add 2014-01-20
	const int lgd = GridT.lgd;

	double** ppx;
	double** ppy;
	double** ppz;

	if(GridT.lgd>0)
	{
		ppx = new double*[lgd];
	  	ppy = new double*[lgd];
	  	ppz = new double*[lgd];
	  	for(int i=0; i<lgd; i++)
	  	{
			ppx[i] = new double[lgd];
			ppy[i] = new double[lgd];
			ppz[i] = new double[lgd];
			ZEROS( ppx[i], lgd);
			ZEROS( ppy[i], lgd);
			ZEROS( ppz[i], lgd);
		}
	}
	
	Vector3<double> tau1, dtau;
	for(int T1=0; T1<ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for(int I1=0; I1< atom1->na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);
			if(GridT.in_this_processor[iat])
			{
				assert( lgd > 0 );

				const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				// get the start positions of elements.
				const int DM_start = LNNR.nlocstartg[iat];
				// get the coordinates of adjacent atoms.
				tau1 = atom1->tau[I1];
				GridD.Find_atom(tau1);
				// search for the adjacent atoms.
				int nad = 0;
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
				{
					// get iat2
					const int T2 = GridD.getType(ad);
					const int I2 = GridD.getNatom(ad);
					const int iat2 = ucell.itia2iat(T2, I2);
					if(GridT.in_this_processor[iat2])
					{
						Atom* atom2 = &ucell.atoms[T2];
						dtau = GridD.getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * ucell.lat0;
						double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
						if(distance < rcut)
						{
							const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
							int ixxx = DM_start + LNNR.find_R2st[iat][nad];
							for(int iw=0; iw<atom1->nw; iw++)
							{
								const int iw_all = start1+iw;
								const int iw_local = GridT.trace_lo[iw_all];
								// iw1_lo
								double *vijx = ppx[iw_local];
								double *vijy = ppy[iw_local];
								double *vijz = ppz[iw_local];

								double *vRx = &pvdpx[ixxx]; //just fold R to normal matrix.
								double *vRy = &pvdpy[ixxx];
								double *vRz = &pvdpz[ixxx];

								int* iw2_lo = &GridT.trace_lo[start2];
								int* iw2_end = iw2_lo + atom2->nw;

								for(; iw2_lo<iw2_end; ++iw2_lo, ++vRx, ++vRy, ++vRz)
								{
									vijx[iw2_lo[0]] += vRx[0] ;
									vijy[iw2_lo[0]] += vRy[0] ;
									vijz[iw2_lo[0]] += vRz[0] ;
								}
								ixxx += atom2->nw;
							}
							++nad;
						}//end distance<rcut
					}
				}//end ad
			}
		}//end ia
	}//end it



	double* tmp = new double[NLOCAL*3];
	for(int i=0; i<NLOCAL; ++i)
	{
		ZEROS(tmp, 3*NLOCAL);
		const int mug = GridT.trace_lo[i];
		// if the row element is on this processor
		if(mug>=0)
		{
			//ofs_running << " i=" << i << " mug=" << mug << endl;
			for(int j=0; j<NLOCAL; ++j)
			{
				const int nug = GridT.trace_lo[j];
				// if the col element is on this processor
				if(nug>=0)
				{
	//				if(mug<nug)
	//				{
						const int index = 3*j;
						tmp[index] = ppx[mug][nug];
						tmp[index+1] = ppy[mug][nug];
						tmp[index+2] = ppz[mug][nug];
	//				}
	//				else
	//				{
					//	tmpx[j] = 0.0;
					//	tmpy[j] = 0.0;
					//	tmpz[j] = 0.0;
	//				}
				}
			}
		}
		// collect the matrix after folding.
		Parallel_Reduce::reduce_double_pool( tmp, NLOCAL*3 );
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			const int iat = ucell.iwt2iat[i];
			const int index = 3*j;
			fvl_dphi[iat][0] += 2.0*tmp[index];	
			fvl_dphi[iat][1] += 2.0*tmp[index+1];	
			fvl_dphi[iat][2] += 2.0*tmp[index+2];	
		}
	}
	delete[] tmp;

	// mohan add 2014-01-20
	if(GridT.lgd > 0)
	{
		//-------------------------
		// delete the tmp matrix.
		//-------------------------
		for(int i=0; i<GridT.lgd; i++)
		{
			delete[] ppx[i];
			delete[] ppy[i];
			delete[] ppz[i];
		}
		delete[] ppx;
		delete[] ppy;
		delete[] ppz;
	}
	timer::tick("Gint_k","folding_force");
	return;
}



// folding the matrix for 'ik' k-point.
// H(k)=\sum{R} H(R)exp(ikR) 
void Gint_k::folding_vl_k(const int &ik)
{
	TITLE("Gint_k","folding_vl_k");
	timer::tick("Gint_k","folding_vl_k",'G');

	if(!pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
	}

	//####################### EXPLAIN #################################
	// 1. what is GridT.lgd ?
	// GridT.lgd is the number of orbitals in each processor according
	// to the division of real space FFT grid.
	// 
	// 2. why the folding of vlocal is different from folding of 
	// < phi_0i | T+Vnl | phi_Rj > ?
	// Because the (i,j) is different for T+Vna and Vlocal+Vna
	// The first part is due to 2D division of H and S matrix,
	// The second part is due to real space division. 
	// 
	// here we construct a temporary matrix to store the
	// matrix element < phi_0 | Vlocal+Vna | phi_R >
	// Vna appears when we use it.
	//#################################################################
	this->ik_now = ik;
	this->pvp = new complex<double>*[GridT.lgd];
	for(int i=0; i<GridT.lgd; i++)
	{
		this->pvp[i] = new complex<double>[GridT.lgd];
		ZEROS( this->pvp[i], GridT.lgd);
	}

	if(!reduced)
	{	
		Vector3<double> dR;
		double arg;
		complex<double> phase;
		complex<double> *pp1;
		double *pp2;
		int count;
		for(int k=0; k<GridT.nutot; k++)
		{
			const int R1x = GridT.ucell_index2x[k];
			const int R1y = GridT.ucell_index2y[k];
			const int R1z = GridT.ucell_index2z[k];

			const int dimk = GridT.lgd*k;
			for(int m=0; m<GridT.nutot; m++)
			{
				//------------------------------------------------
				// exp(k dot dR)
				// dR is the index of box in Crystal coordinates
				//------------------------------------------------
				dR.x = GridT.ucell_index2x[m] - R1x;
				dR.y = GridT.ucell_index2y[m] - R1y;
				dR.z = GridT.ucell_index2z[m] - R1z;

				arg = (kv.kvec_d[ this->ik_now ] * dR) * TWO_PI;
				phase = complex<double>(cos(arg), sin(arg));
				for(int i=0; i<GridT.lgd; i++)
				{
					pp1 = this->pvp[i];
					pp2 = this->pvpR[i+GridT.lgd*m];
					count = dimk;
					for(int j=0; j<GridT.lgd; j++)
					{
						// folding matrix
						pp1[j] += pp2[count] * phase;
						++count;
					}
				}
			}
		}
	}
	else
	{
		int lgd = 0;
		Vector3<double> tau1, dtau, dR;
		for(int T1=0; T1<ucell.ntype; ++T1)
		{
			for(int I1=0; I1<ucell.atoms[T1].na; ++I1)
			{
				// get iat
				const int iat = ucell.itia2iat(T1,I1);
				// atom in this grid piece.
				if(GridT.in_this_processor[iat])
				{
					Atom* atom1 = &ucell.atoms[T1];
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

					// get the start positions of elements.
					const int DM_start = LNNR.nlocstartg[iat];

					// get the coordinates of adjacent atoms.
					tau1 = ucell.atoms[T1].tau[I1];
					GridD.Find_atom(tau1);	
					// search for the adjacent atoms.
					int nad = 0;


					for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
					{
						// get iat2
						const int T2 = GridD.getType(ad);
						const int I2 = GridD.getNatom(ad);
						const int iat2 = ucell.itia2iat(T2, I2);


						// adjacent atom is also on the grid.
						if(GridT.in_this_processor[iat2])
						{
							Atom* atom2 = &ucell.atoms[T2];
							dtau = GridD.getAdjacentTau(ad) - tau1;
							double distance = dtau.norm() * ucell.lat0;
							double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

							// for the local part, only need to calculate <phi_i | phi_j> within range
							// mohan note 2012-07-06
							if(distance < rcut)
							{
								const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

								// calculate the distance between iat1 and iat2.
								// Vector3<double> dR = GridD.getAdjacentTau(ad) - tau1;
								dR.x = GridD.getBox(ad).x;
								dR.y = GridD.getBox(ad).y;
								dR.z = GridD.getBox(ad).z;

								// calculate the phase factor exp(ikR).
								const double arg = (kv.kvec_d[ this->ik_now ] * dR) * TWO_PI;
								complex<double> phase = complex<double>(cos(arg), sin(arg));
								int ixxx = DM_start + LNNR.find_R2st[iat][nad];
								for(int iw=0; iw<atom1->nw; iw++)
								{
									// iw1_lo
									complex<double> *vij = this->pvp[GridT.trace_lo[start1+iw]];


									int* iw2_lo = &GridT.trace_lo[start2];
									int* iw2_end = iw2_lo + atom2->nw;

									if(VNA)
									{
										// get the <phi | V | phi>(R) Hamiltonian.
										double *vijR = &pvpR_reduced[ixxx];
										double *vijR2 = &pvnapR_reduced[ixxx];
										for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR, ++vijR2)
										{
											vij[iw2_lo[0]] += ( vijR[0] + vijR2[0] ) * phase; 
										}
									}
									else
									{
										// get the <phi | V | phi>(R) Hamiltonian.
										double *vijR = &pvpR_reduced[ixxx];
										for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR)
										{
											vij[iw2_lo[0]] += vijR[0] * phase; 
										}
									}
									ixxx += atom2->nw;
									++lgd;
								}

								++nad;
							}// end distane<rcut
						}
					}// end ad
				}
			}// end ia
		}// end it

		//------------------
		// To test the pvpR
		//------------------
/*
		for(int i=0; i<LNNR.nlocdimg[0]; i++)
		{
			const int DM_start = LNNR.nlocstartg[0];
			const int j = i + DM_start;
			if( abs(pvpR_reduced[j]) > 1.0e-5  )
			{
//				cout << " pvpR_reduced[" << i <<"] = " << pvpR_reduced[j] << endl;
			}
		}
*/

	}

	//----------------------
	// Print the pvp matrix
	//----------------------
/*
	cout << " pvp matrix:" << endl;
	for(int i=0; i<GridT.lgd; i++)
	{
		for(int j=0; j<GridT.lgd; j++)
		{
			cout << setw(15) << pvp[i][j].real();
		}
		cout << endl;
	}
	*/

	// Distribution of data.
	timer::tick("Gint_k","Distri",'G');
	complex<double>* tmp;
	for (int i=0; i<NLOCAL; i++)
	{
		tmp = new complex<double>[NLOCAL];
		ZEROS(tmp, NLOCAL);
		const int mug = GridT.trace_lo[i];
		// if the row element is on this processor.
		if (mug >= 0)
		{
			for (int j=0; j<NLOCAL; j++)
			{
				const int nug = GridT.trace_lo[j];
				// if the col element is on this processor.
				if (nug >=0)
				{
					if (mug <= nug)
					{
						// pvp is symmetric, only half is calculated.
						tmp[j] = this->pvp[mug][nug];
					}
					else
					{
						// need to get elements from the other half.
						// I have question on this! 2011-02-22
						tmp[j] = conj(this->pvp[nug][mug]);
					}
				}
			}
		}
		// collect the matrix after folding.
		Parallel_Reduce::reduce_complex_double_pool( tmp, NLOCAL );

		//-----------------------------------------------------
		// NOW! Redistribute the Hamiltonian matrix elements
		// according to the HPSEPS's 2D distribution methods.
		//-----------------------------------------------------
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			// set the matrix value.
			LM.set_HSk(i,j,tmp[j],'L');
		}
		delete[] tmp;
	}

	// delete the tmp matrix.
	for(int i=0; i<GridT.lgd; i++)
	{
		delete[] pvp[i];
	}
	delete[] pvp;
	timer::tick("Gint_k","Distri",'G');

	timer::tick("Gint_k","folding_vl_k",'G');
	return;
}

void Gint_k::set_ijk_atom(const int &grid_index, const int &size,
	double*** psir_ylm, double*** dr, bool** cal_flag, 
	double** distance, double* ylma, const double &delta_r)
{
	const Numerical_Orbital_Lm* pointer;
	double mt[3];
	for (int id=0; id<size; id++)
	{
		// (2.1) get the atom type and atom index.
		const int mcell_index = GridT.bcell_start[grid_index] + id;	
		const int imcell = GridT.which_bigcell[mcell_index];
		const int iat = GridT.which_atom[mcell_index];
		const int it = ucell.iat2it[ iat ];
		const int ia = ucell.iat2ia[ iat ];

		// (2.2) get the distance between the grid and the atom.
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
			+ dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);

			//if(distance[ib][id] <= ORB.Phi[it].getRcut())
			// mohan reset this 2012-06-27
			// if(distance[ib][id] < ORB.Phi[it].getRcut() )
			// mohan reset this 2013-07-02 in Princeton
			// we should make absolutely sure that the distance is smaller than ORB.Phi[it].getRcut
			// this should be consistant with LCAO_nnr::cal_nnrg function 
			// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
			// distance = 7.0000000000000000
			// ORB.Phi[it].getRcut = 7.0000000000000008
			if(distance[ib][id] < (ORB.Phi[it].getRcut() - 1.0e-15) )
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

			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position = distance[ib][id] / delta_r;

			const int ip = static_cast<int>(position);

			const double dx = position - ip;
			const double dx2 = dx * dx;
			const double dx3 = dx2 * dx;

			const double c3 = 3.0*dx2-2.0*dx3;
			const double c1 = 1.0-c3;
			const double c2 = (dx-2.0*dx2+dx3)*delta_r;
			const double c4 = (dx3-dx2)*delta_r;

			Atom* atom1 = &ucell.atoms[it];
			double tmp=0.0;//mohan fix bug 2011-05-04
			for (int iw=0; iw< atom1->nw; iw++)
			{
				if ( atom1->iw2_new[iw] )
				{
					pointer = &ORB.Phi[it].PhiLN(
							atom1->iw2l[iw],
							atom1->iw2n[iw]);

					// Efficient!! to get the orbital value at this point.
					tmp = c1*pointer->psi_uniform[ip] 
						+ c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] 
						+ c4*pointer->dpsi_uniform[ip+1];
				}
				psir_ylm[ib][id][iw] = tmp * ylma[atom1->iw2_ylm[iw]];
			}// end iw.
		}//end ib
	}// int id
	return;
}


void Gint_k::set_ijk_atom_vna(const int &grid_index, const int &size,
	double*** psir_ylm, double*** dr, bool** cal_flag, 
	double** distance, double* ylma, const double &delta_r,
	const Grid_Technique &gt, double* vna3d)
{
	const Numerical_Orbital_Lm* pointer;
	double mt[3];
	for (int id=0; id<size; id++)
	{
		// (2.1) get the atom type and atom index.
		const int mcell_index = gt.bcell_start[grid_index] + id;	
		const int imcell = gt.which_bigcell[mcell_index];
		const int iat = gt.which_atom[mcell_index];
		const int it = ucell.iat2it[ iat ];
		const int ia = ucell.iat2ia[ iat ];

		// (2.2) get the distance between the grid and the atom.
		mt[0] = gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0];
		mt[1] = gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1];
		mt[2] = gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2];

		// as big as 64!
		for(int ib=0; ib<gt.bxyz; ib++)
		{
			// meshcell_pos: z is the fastest
			dr[ib][id][0] = gt.meshcell_pos[ib][0] + mt[0];
			dr[ib][id][1] = gt.meshcell_pos[ib][1] + mt[1];
			dr[ib][id][2] = gt.meshcell_pos[ib][2] + mt[2];

			distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
			+ dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);


			//mohan fix bug 2012-07-02, should be >=
			if(distance[ib][id] >= ORB.Vna[it].rcut &&
		   	   distance[ib][id] >= ORB.Phi[it].getRcut())
			{
				cal_flag[ib][id]=false;
//				ofs_running << " ibbb=" << ib << " iddd=" << id << " cal_flag=" << cal_flag[ib][id] << endl;
				continue;
			}

			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;

			

			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position = distance[ib][id] / delta_r;

			const int iq = static_cast<int>(position);
			const double x0 = position - static_cast<double>(iq);
			const double x1 = 1.0 - x0;
			const double x2 = 2.0 - x0;
			const double x3 = 3.0 - x0;
			const double x12 = x1*x2/6.0;
			const double x03 = x0*x3/2.0;



			// mohan add 2012-06-13
			double ccc;
			if(distance[ib][id] <= ORB.Vna[it].rcut)
			{
				ccc = (x12*(ORB.Vna[it].vna_u[iq]*x3
				+ORB.Vna[it].vna_u[iq+3]*x0)
				+ x03*(ORB.Vna[it].vna_u[iq+1]*x2
				-ORB.Vna[it].vna_u[iq+2]*x1));
				vna3d[ib] += ccc;
			}


			//if(distance[ib][id] <= ORB.Phi[it].getRcut())
			// mohan fix bug 2012-06-27
			// if use '<=', there may be atoms out of plan.
			if(distance[ib][id] < ORB.Phi[it].getRcut())
			{
				cal_flag[ib][id]=true;
//				ofs_running << " ibbb=" << ib << " iddd=" << id << " cal_flag=" << cal_flag[ib][id] << endl;
			}
			else
			{
				cal_flag[ib][id]=false;
//				ofs_running << " ibbb=" << ib << " iddd=" << id << " cal_flag=" << cal_flag[ib][id] << endl;
				continue;
			}



			Ylm::sph_harm ( ucell.atoms[it].nwl,
					dr[ib][id][0] / distance[ib][id],
					dr[ib][id][1] / distance[ib][id],
					dr[ib][id][2] / distance[ib][id],
					ylma);





			const int ip = static_cast<int>(position);

			const double dx = position - ip;
			const double dx2 = dx * dx;
			const double dx3 = dx2 * dx;

			const double c3 = 3.0*dx2-2.0*dx3;
			const double c1 = 1.0-c3;
			const double c2 = (dx-2.0*dx2+dx3)*delta_r;
			const double c4 = (dx3-dx2)*delta_r;







			Atom* atom1 = &ucell.atoms[it];
			double tmp=0.0;//mohan fix bug 2011-05-04
			for (int iw=0; iw< atom1->nw; iw++)
			{
				if ( atom1->iw2_new[iw] )
				{
					pointer = &ORB.Phi[it].PhiLN(
							atom1->iw2l[iw],
							atom1->iw2n[iw]);

					// Efficient!! to get the orbital value at this point.
					tmp = c1*pointer->psi_uniform[ip] 
						+ c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] 
						+ c4*pointer->dpsi_uniform[ip+1];
				}
				psir_ylm[ib][id][iw] = tmp * ylma[atom1->iw2_ylm[iw]];
//				psir_ylm[ib][id][iw] = 1.0; 
			}// end iw.
		}//end ib
	}// int id



	/*
	if(grid_index==0)
	{
		for(int id=0; id<size; ++id)
		{
        const int mcell_index = GridT.bcell_start[grid_index] + id;
        const int imcell = GridT.which_bigcell[mcell_index];
        const int iat = GridT.which_atom[mcell_index];
        const int it = ucell.iat2it[ iat ];

			for(int ib=0; ib<pw.bxyz; ++ib)
			{
				ofs_running << " atom=" << id << " ib=" << ib << " it=" << ucell.atoms[it].label 
				<< " distance=" << distance[ib][id] 
				<< " cal_flag=" << cal_flag[ib][id] << endl;	
			}
		} 
	}
	*/

	return;
}



void Gint_k::set_ijk_atom_fvna(const int &grid_index, const int &size,
    double*** psir_ylm, double*** dr, bool** cal_flag,
    double** distance, double* ylma, const double &delta_r,
    double*** dphi_x, double*** dphi_y, double*** dphi_z,
	const Grid_Technique &gt, double* vna3d)
{
    const Numerical_Orbital_Lm* pointer;
    double mt[3];
    // Peize Lin change rly, grly 2016-08-26
    static vector<double> rly;
    static vector<vector<double>> grly;
    for (int id=0; id<size; id++)
    {
        // (2.1) get the atom type and atom index.
        const int mcell_index = gt.bcell_start[grid_index] + id;
        const int imcell = gt.which_bigcell[mcell_index];
        const int iat = gt.which_atom[mcell_index];
        const int it = ucell.iat2it[ iat ];
        const int ia = ucell.iat2ia[ iat ];
        Atom *atom = &ucell.atoms[it];

        // (2.2) get the distance between the grid and the atom.
        mt[0] = gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0];
        mt[1] = gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1];
        mt[2] = gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2];

        for(int ib=0; ib<gt.bxyz; ib++) //mohan fix bug 2012-06-29 
        {
            // meshcell_pos: z is the fastest
            dr[ib][id][0] = gt.meshcell_pos[ib][0] + mt[0];
            dr[ib][id][1] = gt.meshcell_pos[ib][1] + mt[1];
            dr[ib][id][2] = gt.meshcell_pos[ib][2] + mt[2];

            distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0]
            + dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);


			// mohan fix bug 2012-07-02
            if(distance[ib][id] >= ORB.Vna[it].rcut &&
               distance[ib][id] >= ORB.Phi[it].getRcut())
            {
                cal_flag[ib][id]=false;
                continue;
            }

// don't need, deal with in the following
//			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;


            // get the 'phi' and 'dphi'.
            Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[ib][id][0], dr[ib][id][1], dr[ib][id][2], rly, grly);

//          if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
//          Ylm::sph_harm ( ucell.atoms[it].nwl,
//                  dr[ib][id][0] / distance[ib][id],
//                  dr[ib][id][1] / distance[ib][id],
//                  dr[ib][id][2] / distance[ib][id],
//                  ylma);

            // these parameters are about interpolation
            // because once we know the distance from atom to grid point,
            // we can get the parameters we need to do interpolation and
            // store them first!! these can save a lot of effort.
            const double position = distance[ib][id] / delta_r;

            const int iq = static_cast<int>(position);
            const double x0 = position - static_cast<double>(iq);
            const double x1 = 1.0 - x0;
            const double x2 = 2.0 - x0;
            const double x3 = 3.0 - x0;
            const double x12 = x1*x2 / 6.0;
            const double x03 = x0*x3 / 2.0;


            // mohan add 2012-06-13
            double ccc;
            if(distance[ib][id] <= ORB.Vna[it].rcut)
            {
                ccc = (x12*(ORB.Vna[it].vna_u[iq]*x3
                +ORB.Vna[it].vna_u[iq+3]*x0)
                + x03*(ORB.Vna[it].vna_u[iq+1]*x2
                -ORB.Vna[it].vna_u[iq+2]*x1));
                vna3d[ib] += ccc;
            }

			
            //if(distance[ib][id] <= ORB.Phi[it].getRcut())
			// mohan fix 2012-06-27
            if(distance[ib][id] < ORB.Phi[it].getRcut())
            {
                cal_flag[ib][id]=true;
            }
            else
            {
                cal_flag[ib][id]=false;
                continue;
            }			


			double tmp = 0.0;//mohan fix bug 2011-05-04
            double dtmp = 0.0;
            for (int iw=0; iw< atom->nw; iw++)
            {
                if ( atom->iw2_new[iw] )
                {
                    pointer = &ORB.Phi[it].PhiLN(
                            atom->iw2l[iw],
                            atom->iw2n[iw]);

                    if(iq >= pointer->nr_uniform-4)
                    {
                        tmp = dtmp = 0.0;
                    }
                    else
                    {
                        // Efficient!! to get the orbital value at this point.
                        tmp = x12*(pointer->psi_uniform[iq]*x3
                            +pointer->psi_uniform[iq+3]*x0)
                        + x03*(pointer->psi_uniform[iq+1]*x2
                        -pointer->psi_uniform[iq+2]*x1);

                        dtmp = x12*(pointer->dpsi_uniform[iq]*x3
                        +pointer->dpsi_uniform[iq+3]*x0)
                        + x03*(pointer->dpsi_uniform[iq+1]*x2
                        -pointer->dpsi_uniform[iq+2]*x1);

                        //dtmp = x12*(pointer->psi_uniform[iq]*x3
                        //+pointer->psi_uniform[iq+3]*x0)
                        //+ x03*(pointer->psi_uniform[iq+1]*x2
                        //-pointer->psi_uniform[iq+2]*x1);
                    }
                }// new l is used.

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
                        dphi_x[ib][id][iw] = dphi_y[ib][id][iw] = dphi_z[ib][id][iw] = 0.0;
                    }
                    else
                    {
                        pointer = &ORB.Phi[it].
                            PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

                        double Zty = pointer->zty;
                        psir_ylm[ib][id][iw] = Zty * rly[idx_lm];
                        dphi_x[ib][id][iw] = Zty * grly[idx_lm][0];
                        dphi_y[ib][id][iw] = Zty * grly[idx_lm][1];
                        dphi_z[ib][id][iw] = Zty * grly[idx_lm][2];
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
                    dphi_x[ib][id][iw] = tmpdphi_rly * dr[ib][id][0]  + tmprl * grly[idx_lm][0];
                    dphi_y[ib][id][iw] = tmpdphi_rly * dr[ib][id][1]  + tmprl * grly[idx_lm][1];
                    dphi_z[ib][id][iw] = tmpdphi_rly * dr[ib][id][2]  + tmprl * grly[idx_lm][2];
                }
            }
        }//end ib
    }// int id
    return;
}




