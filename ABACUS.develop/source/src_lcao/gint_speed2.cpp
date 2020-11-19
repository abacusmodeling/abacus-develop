#include "gint_speed.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"

void Gint_Speed::save_phi(void)
{
    TITLE("Gint_Speed","save_phi");
	timer::tick("Gint_Speed","save_phi",'k');
	cout << " SAVE ORBITALS ON GRID." << endl;

	const int nat = ucell.nat; 
	ofs_running << " allocate_phiylm = " << allocate_phiylm << endl;

	if(allocate_phiylm)
	{
		delete[] vindex;
		for(int iat=0; iat<nat; ++iat)
		{
			if( nr[iat] > 0)
			{
				delete[] grid_label[iat];
			}
		}
		delete[] grid_label;

		for(int iat=0; iat<nat; ++iat)
		{
			if( nr[iat] > 0 )
			{
				const int it = ucell.iat2it[iat];
				const int nw = ucell.atoms[it].nw;
				for(int iw=0; iw<nw; ++iw)
				{
					delete[] phiylm[iat][iw];
				}
				delete[] phiylm[iat];
			}
		}
		delete[] phiylm;

		
		// nr delete at last!
		delete[] nr;
	}


	const int nbx = GridT.nbx;
	const int nby = GridT.nby;
	const int nbz_start = GridT.nbzp_start;
	const int nbz = GridT.nbzp;

	const int bx = pw.bx;
	const int by = pw.by;
	const int bz = pw.bz;
	const int bxyz = pw.bxyz;

	const int ncx = pw.ncx;
	const int ncy = pw.ncy;
	const int ncz = pw.ncz;
	const int nczp = pw.nczp;
	const int nczp_start = pw.nczp_start;

//	cout << " nbx=" << nbx << " nby=" << nby << " nbz=" << nbz << " nbz_start=" << nbz_start << endl;
//	cout << " bx=" << bx << " by=" << by << " bz=" << bz << " nbz=" << nbz << endl;
//	cout << " ncx=" << ncx << " ncy=" << ncy << " ncz=" << ncz << " nczp=" << nczp << " nczp_start=" << nczp_start << endl;

	this->nr = new int[nat];
	int *nr2 = new int[nat];
	ZEROS(nr, nat);
	ZEROS(nr2, nat);

	//-----------------------------------
	// PART1: get the nr for each atom.
	//-----------------------------------
	for (int k=nbz_start; k<nbz_start+nbz; ++k) // FFT grid
	{
		for (int j=0; j<nby; ++j)
		{
			for (int i=0; i< nbx; ++i)
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
					double mt0 = GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
					double mt1 = GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
					double mt2 = GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];
					for(int ib=0; ib<bxyz; ++ib)
					{
						double dr0 = GridT.meshcell_pos[ib][0] + mt0;
						double dr1 = GridT.meshcell_pos[ib][1] + mt1;
						double dr2 = GridT.meshcell_pos[ib][2] + mt2;
						double distance = std::sqrt(dr0*dr0+dr1*dr1+dr2*dr2);
						if(distance <= ORB.Phi[it].getRcut())
						{
							++nr[iat];
						}
					}
				}
			}
		}
	}


	




	//--------------------------------------------------
	// PART2: record the wave functions for each atom.
	//--------------------------------------------------
	// I, nw, r 
	int cm = 0;
	this->phiylm = new double**[nat];
	for(int iat=0; iat<nat; ++iat)
	{
		const int nr_now = nr[iat];
		if(nr_now > 0) // need to be updated later, mohan 2012-06-27
		{
			const int it = ucell.iat2it[iat];
			phiylm[iat] = new double*[ucell.atoms[it].nw];
			
	//		ofs_running << " iat=" << iat << " nr_now=" << nr_now << endl;
			for(int iw=0; iw<ucell.atoms[it].nw; ++iw)
			{
				phiylm[iat][iw] = new double[nr_now];
				ZEROS( phiylm[iat][iw], nr_now);
				cm += nr_now;
			}
		}
	}


	this->grid_label = new int*[nat];
	
	for(int iat=0; iat<nat; ++iat)
	{
		const int nr_now = this->nr[iat];
		if(nr_now > 0)
		{
			grid_label[iat] = new int[nr[iat]];
			ZEROS(grid_label[iat], nr[iat]);
		}
	}

	this->vindex = new int[pw.nrxx];
	ZEROS(vindex, pw.nrxx);


	this->allocate_phiylm = true;

//	ofs_running << " cm=" << cm << endl;
//	ofs_running << " bxyz=" << bxyz << endl;
	ofs_running << " Phi Memory is " <<  Memory::record("Grid_Int", "phiylm", cm*bxyz, "double") << " MB" << endl; 	

	// get the wave functions on the x-y plane
	int nnnmax=0;
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax = max(nnnmax, nnn[T]);
	}

	const Numerical_Orbital_Lm* pointer;
	const double delta_r = ORB.dr_uniform;
	double* ylma = new double[nnnmax]; // Ylm for each atom: [bxyz, nnnmax]
	double mt[3]={0,0,0};	

	int ir_bigbox = 0;
	for (int k=nbz_start; k<nbz_start+nbz; k++) // FFT grid
	{
		for (int j=0; j<nby; j++)
		{
			for (int i=0; i< nbx; i++)
			{
				this->grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
				const int size = GridT.how_many_atoms[ this->grid_index ];
				
				if(size==0) continue;

				int ir_now = ir_bigbox; 
				for(int ib=0; ib<bxyz; ib++)
				{
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

						double dr0 = GridT.meshcell_pos[ib][0] + mt[0];
						double dr1 = GridT.meshcell_pos[ib][1] + mt[1];
						double dr2 = GridT.meshcell_pos[ib][2] + mt[2];
						
						double distance = std::sqrt(dr0*dr0+dr1*dr1+dr2*dr2);
						if(distance <= ORB.Phi[it].getRcut())
						{
						}
						else
						{
							continue;
						}
						if (distance < 1.0E-9) distance += 1.0E-9;
						
						Ylm::sph_harm ( ucell.atoms[it].nwl,
								dr0 / distance,
								dr1 / distance,
								dr2 / distance,
								ylma);
						const double position = distance / delta_r;
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
						for (int iw=0; iw< atom1->nw; iw++)
						{
							double phi;
							if ( atom1->iw2_new[iw] )
							{
								pointer = &ORB.Phi[it].PhiLN(
										atom1->iw2l[iw],
										atom1->iw2n[iw]);
								phi = c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
									+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
							}
							phiylm[iat][iw][nr2[iat]] = phi * ylma[atom1->iw2_ylm[iw]];
						}

						// recall the global index of each grid for each orbital
						// ir_bigbox = global index(i, j, k)

						//-----------------------------------------------------
						// (1)
						// Be careful about one special case!
						// if grid_label[iat][nr_now] == grid_label[iat][nr_now+1],
						// means, there are at least two 'identical' atoms
						// have non-zero wave functios on this grid. 
						// (2)
						// Because we are now circleing for ir_bigbox, so if 
						// there appears same 'iat', 
						// nr2[iat] is blind between 'different atom (same iat) on same grid'
						// and 'same atom on different grid'
						//-----------------------------------------------------
						//assert(nr[iat]>0);
						grid_label[iat][nr2[iat]] = ir_now;
						++nr2[iat];
					}//end size
					++ir_now;
				}//end ibxyz


				ir_now = ir_bigbox;
				for(int ii=0; ii<bx; ii++)
				{
					const int iii = i*bx + ii;
					const int ipart = iii*ncy*nczp;
					for(int jj=0; jj<by; jj++)
					{
						const int jjj = j*by + jj;
						const int jpart = jjj*nczp;
						for(int kk=0; kk<bz; kk++)
						{
							const int kkk = k*bz + kk;

							//assert(ir_now<pw.nrxx);
							vindex[ir_now] = (kkk-nczp_start) + jpart + ipart;
							++ir_now;
							//assert(vindex[bindex] < pw.nrxx);
						}
					}
				}

				ir_bigbox += bxyz;
			}//end i
		}//end j
	}//end k

	delete[] ylma;



	/*
	ofs_running << setw(10) << "index" << setw(25) << "nr" << setw(25) << "nr2" << endl;
	for(int iat=0; iat<nat; ++iat)
	{
		assert(nr[iat]==nr2[iat]);
		ofs_running << setw(10) << iat+1 << setw(25) << nr[iat] << setw(25) << nr2[iat] << endl;
	}
	*/
	
	//----------------------------------------------------
	// PART3: calculate the dimension of tmporary arrays
	//----------------------------------------------------

	//----------------------------------------------------
	// nov is the maximal overlap mesh number between
	// atom 1 and atom 2.
	//----------------------------------------------------
	this->nov=0;
	for(int iat=0; iat<nat; ++iat)
	{
		const int nr_now = this->nr[iat];
		if(nr_now == 0) continue;

		const int it = ucell.iat2it[ iat ];
		const int ia = ucell.iat2ia[ iat ];
		const int start1 = ucell.itiaiw2iwt(it, ia, 0);
		const int iw_lo = GridT.trace_lo[start1];
		if(iw_lo<0)continue;
		for(int iat2=iat; iat2<nat; ++iat2)
		{
			const int nr2_now = this->nr[iat2];
			if(nr2_now == 0) continue;

			const int it2 = ucell.iat2it[ iat2 ];
			const int ia2 = ucell.iat2ia[ iat2 ];
			const int start2 = ucell.itiaiw2iwt(it2, ia2, 0);
			const int iw2_lo = GridT.trace_lo[start2];
			if(iw2_lo<0)continue;



			// for example
			// the positive number is the grid point.
			// positive number is always increasing.
			// psi1: 1 3 3 5 7
			// psi2: 2 3 3 5 9
			int is=0;
			int is_rec=0;
			int ir2_start=0;//start ir for atom2
			int ee=0; //count the overlap grid number.
			bool same=false;


			// searching in the range of atom 1.
			for(int ir=0; ir< nr_now; ++ir)
			{
				if( ir>0 )
				{
					//---------------------------------------
					// means 'idential' atoms have at least
					// twice phi on this grid.
					//---------------------------------------
					if( grid_label[iat][ir] == grid_label[iat][ir-1] )
					{
						// it's the same grid.
						same = true;
						// again, searching begins from the previous one.
						// for example, the second 3 in atom 1,
						// must search from the previous one. 
						ir2_start = is_rec;	
					}
					else
					{
						// start from the bigger one, when
						// grid_label[iat2][j] > grid_label[iat][i] 
						same = false;
						ir2_start = is;
					}
				}

				// is_rec: 
				// stop to begin record how many 'identical'
				// atoms have Phi on this grid.
				is_rec = ir2_start;

				//------------------------------------------
				// (1) compare betwen different atoms!
				// (2) ee++ when iat and iat2 are on the same grid.
				// (3) can't miss any grid_label
				//------------------------------------------
				for(int jr=ir2_start; jr<nr2_now; ++jr)
				{
					// this means another grid appears,
					// so we can break the search,
					// 
					// the relationship between grid_label[iat2][j]
					// and grid_label[iat][i+1] is stil no known,
					// so we record the position j (belongs to atom2). 
					// the recorded index is in 'is' variable.

					if( grid_label[iat2][jr] > grid_label[iat][ir] )
					{
						is=jr;
						break;
					}
					else if( grid_label[iat2][jr] == grid_label[iat][ir])
					{
						++ee;
					}
				}// j
			}// i
//			ofs_running << "iat=" << iat << " iat2=" << iat2 << " ee=" << ee << endl;
			this->nov = std::max(nov, ee);
		}//iat2
	}//iat


	ofs_running << " nov=" << nov << " nrxx=" << pw.nrxx << endl;
//	cout << " nov=" << nov << " nrxx=" << pw.nrxx << endl;

#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	timer::tick("Gint_Speed","save_phi",'k');
    return;
}

