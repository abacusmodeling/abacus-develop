#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "blas_interface.h"

//#include <vector>

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
}

inline void setVindex(const int ncyz, const int ibx, const int jby, const int kbz, 
					int* vindex)
{				
	int bindex=0;
	// z is the fastest, 
	for(int ii=0; ii<pw.bx; ii++)
	{
		const int ipart=(ibx + ii) * ncyz + kbz;
		for(int jj=0; jj<pw.by; jj++)
		{
			const int jpart=(jby + jj) * pw.nczp + ipart;
			for(int kk=0; kk<pw.bz; kk++)
			{
				vindex[bindex]=kk + jpart; 
				++bindex;
			}
		}
	}
}

inline void cal_psir_ylm(int size, int grid_index, double delta_r,
						double** distance, double* ylma,
						int* at, int* block_index, int* block_iw, int* block_size, 
						int** cal_flag, double** psir_ylm)
{
	const Numerical_Orbital_Lm *pointer;
	double mt[3];
	double dr[3];
	block_index[0]=0;
	for (int id=0; id<size; id++)
	{
		// there are two parameters we want to know here:
		// in which bigcell of the meshball the atom in?
		// what's the cartesian coordinate of the bigcell?
		const int mcell_index=GridT.bcell_start[grid_index] + id;
		const int imcell=GridT.which_bigcell[mcell_index];

		const int iat=GridT.which_atom[mcell_index];
		at[id]=iat;
		
		const int it=ucell.iat2it[iat];
		const int ia=ucell.iat2ia[iat];
		const int start=ucell.itiaiw2iwt(it, ia, 0);
		block_iw[id]=GridT.trace_lo[start];
		Atom* atom=&ucell.atoms[it];
		block_size[id]=atom->nw;
		block_index[id+1]=block_index[id]+atom->nw;
		// meshball_positions should be the bigcell position in meshball
		// to the center of meshball.
		// calculated in cartesian coordinates
		// the vector from the grid which is now being operated to the atom position.
		// in meshball language, is the vector from imcell to the center cel, plus
		// tau_in_bigcell.
		mt[0]=GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
		mt[1]=GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
		mt[2]=GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

		for(int ib=0; ib<pw.bxyz; ib++)
		{
			double *p=&psir_ylm[ib][block_index[id]];
			// meshcell_pos: z is the fastest
			dr[0]=GridT.meshcell_pos[ib][0] + mt[0]; 
			dr[1]=GridT.meshcell_pos[ib][1] + mt[1]; 
			dr[2]=GridT.meshcell_pos[ib][2] + mt[2]; 	

			distance[ib][id]=std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if(distance[ib][id] > (ORB.Phi[it].getRcut()- 1.0e-15)) 
			{
				cal_flag[ib][id]=0;
				ZEROS(p, block_size[id]);
				continue;
			}

			cal_flag[ib][id]=1;
			
			//if(distance[id] > GridT.orbital_rmax) continue;
			//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
			
			Ylm::sph_harm (	ucell.atoms[it].nwl,
					dr[0] / distance[ib][id],
					dr[1] / distance[ib][id],
					dr[2] / distance[ib][id],
					ylma);
			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position=distance[ib][id] / delta_r;
			const int ip=static_cast<int>(position);
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
					pointer=&ORB.Phi[it].PhiLN(
							atom->iw2l[iw],
							atom->iw2n[iw]);
					phi=c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
				}
				*p=phi * ylma[atom->iw2_ylm[iw]];
			} // end iw
		}// end ib
	}// end id
}

inline void cal_meshball_vlocal(int size, int LD_pool, int* block_iw, int* bsize, int* colidx, 
							int** cal_flag, double* vldr3, double** psir_ylm, double** psir_vlbr3, 
							int* vindex, int lgd_now, double** GridVlocal)
{
	char transa='N', transb='T';
	double alpha=1, beta=1;
	
	int allnw=colidx[size];
	for(int ib=0; ib<pw.bxyz; ++ib)
	{
        for(int ia=0; ia<size; ++ia)
        {
            if(cal_flag[ib][ia]>0)
            {
                for(int i=colidx[ia]; i<colidx[ia+1]; ++i)
                {
                    psir_vlbr3[ib][i]=psir_ylm[ib][i]*vldr3[ib];
                }
            }
            else
            {
                for(int i=colidx[ia]; i<colidx[ia+1]; ++i)
                {
                    psir_vlbr3[ib][i]=0;
                }
            }
            
        }
	}

	for(int ia1=0; ia1<size; ++ia1)
	{
		const int iw1_lo=block_iw[ia1];
		int m=bsize[ia1];	
		for(int ia2=0; ia2<size; ++ia2)
		{
			const int iw2_lo=block_iw[ia2];
			if(iw1_lo<=iw2_lo)
			{
                int first_ib=0, last_ib=0;
                for(int ib=0; ib<pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                    {
                        first_ib=ib;
                        break;
                    }
                }
                for(int ib=pw.bxyz-1; ib>=0; --ib)
                {
                    if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                    {
                        last_ib=ib+1;
                        break;
                    }
                }
                int ib_length=last_ib-first_ib;
                if(ib_length<=0) continue;

                int cal_pair_num=0;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    cal_pair_num+=cal_flag[ib][ia1]*cal_flag[ib][ia2];
                }
                
                int n=bsize[ia2];
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &n, &m, &ib_length, &alpha,
                        &psir_vlbr3[first_ib][colidx[ia2]], &LD_pool, 
                        &psir_ylm[first_ib][colidx[ia1]], &LD_pool,  
                        &beta, &GridVlocal[iw1_lo][iw2_lo], &lgd_now);
                }
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                        {
                            int k=1;
                            dgemm_(&transa, &transb, &n, &m, &k, &alpha,
                                &psir_vlbr3[ib][colidx[ia2]], &LD_pool, 
                                &psir_ylm[ib][colidx[ia1]], &LD_pool,  
                                &beta, &GridVlocal[iw1_lo][iw2_lo], &lgd_now);
                        }
                    }
                }
                
			}
		}
	}
}

inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
    //return (localIndex/nblk*nprocs+myproc)*nblk+localIndex%nblk;
}

inline int localIndex(int globalIndex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalIndex%(nblk*nprocs))/nblk);
    return int(globalIndex/(nblk*nprocs))*nblk+globalIndex%nblk;
}

inline int setBufferParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk,
                              int& sender_index_size, int*& sender_local_index, 
                              int*& sender_size_process, int*& sender_displacement_process, 
                              int& sender_size, double*& sender_buffer,
                              int& receiver_index_size, int*& receiver_global_index, 
                              int*& receiver_size_process, int*& receiver_displacement_process, 
                              int& receiver_size, double*& receiver_buffer)
{
    // setup blacs parameters
    int nprows, npcols, nprocs;
    int myprow, mypcol, myproc;
    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    Cblacs_pinfo(&myproc, &nprocs);
    
    // init data arrays
    int current_sender_index_size=GridT.lgd*GridT.lgd*2;
    if(current_sender_index_size > sender_index_size)
    {
        sender_index_size=current_sender_index_size;
        // OUT(ofs_running, "sender_index_size:", sender_index_size);
        delete[] sender_local_index;
        sender_local_index=new int[sender_index_size];
    }

    static bool FIRST_RUN=true;
    if(FIRST_RUN)
    {
        FIRST_RUN=false;
        delete[] sender_size_process;
        sender_size_process=new int[nprocs];
        delete[] sender_displacement_process;
        sender_displacement_process=new int[nprocs];

        delete[] receiver_size_process;
        receiver_size_process=new int[nprocs];
        delete[] receiver_displacement_process;
        receiver_displacement_process=new int[nprocs];
    }

    // build the local index to be sent to other process (sender_local_index),
    //       the global index to be received from other process (receiver_global_index),
    //       the send/receive size/displacement for data exchange by MPI_Alltoall
    int *sender_global_index=new int[current_sender_index_size];

    int pos=0;
    sender_size_process[0]=0;
    for(int iproc=0; iproc<nprocs; ++iproc)
    {
        sender_displacement_process[iproc]=pos;
     
        int iprow, ipcol;
        Cblacs_pcoord(blacs_ctxt, iproc, &iprow, &ipcol);
        
        // find out the global index and local index of elements in each process based on 2D block cyclic distribution
        for(int irow=0, grow=0; grow<NLOCAL; ++irow)
        {
            grow=globalIndex(irow, nblk, nprows, iprow);
            int lrow=GridT.trace_lo[grow];
            if(lrow < 0 || grow >= NLOCAL) continue;
            for(int icol=0, gcol=0; gcol<NLOCAL; ++icol)
            {
                gcol=globalIndex(icol,nblk, npcols, ipcol);
                int lcol=GridT.trace_lo[gcol];
                if(lcol < 0 || gcol >= NLOCAL) continue;
                // if(pos<0 || pos >= current_sender_index_size)
                // {
                //     OUT(ofs_running, "pos error, pos:", pos);
                //     OUT(ofs_running, "irow:", irow);
                //     OUT(ofs_running, "icol:", icol);
                //     OUT(ofs_running, "grow:", grow);
                //     OUT(ofs_running, "gcol:", gcol);
                //     OUT(ofs_running, "lrow:", grow);
                //     OUT(ofs_running, "lcol:", gcol);
                // }
                sender_global_index[pos]=grow;
                sender_global_index[pos+1]=gcol;
                sender_local_index[pos]=lrow;
                sender_local_index[pos+1]=lcol;
                pos+=2;
            }
        }
        sender_size_process[iproc]=pos-sender_displacement_process[iproc];
    }
   
    MPI_Alltoall(sender_size_process, 1, MPI_INT, 
                 receiver_size_process, 1, MPI_INT, comm_2D);

    int current_receiver_index_size=receiver_size_process[0];
    receiver_displacement_process[0]=0;
    for(int i=1; i<nprocs; ++i)
    {
        current_receiver_index_size+=receiver_size_process[i];
        receiver_displacement_process[i]=receiver_displacement_process[i-1]+receiver_size_process[i-1];
    }

    if(current_receiver_index_size > receiver_index_size)
    {
        receiver_index_size=current_receiver_index_size;
        // OUT(ofs_running, "receiver_index_size:", receiver_index_size);
        delete[] receiver_global_index;
        receiver_global_index=new int[receiver_index_size];
    }

    // send the global index in sendBuffer to recvBuffer
    MPI_Alltoallv(sender_global_index, sender_size_process, sender_displacement_process, MPI_INT, 
                  receiver_global_index, receiver_size_process, receiver_displacement_process, MPI_INT, comm_2D);
    
    delete [] sender_global_index;

    // the sender_size_process, sender_displacement_process, receiver_size_process, 
    // and receiver_displacement_process will be used in transfer sender_buffer, which
    // is half size of sender_global_index
    // we have to rebuild the size and displacement for each process
    for (int iproc=0; iproc < nprocs; ++iproc)
    {
        sender_size_process[iproc]=sender_size_process[iproc]/2;
        sender_displacement_process[iproc]=sender_displacement_process[iproc]/2;
        receiver_size_process[iproc]=receiver_size_process[iproc]/2;
        receiver_displacement_process[iproc]=receiver_displacement_process[iproc]/2;
    }
    
    int current_sender_size=current_sender_index_size/2;
    if(current_sender_size > sender_size)
    {
        sender_size=current_sender_size;
        // OUT(ofs_running, "sender_size:", sender_size);
        delete[] sender_buffer;
        sender_buffer=new double[sender_size];
    }
    int current_receiver_size=current_receiver_index_size/2;
    if(current_receiver_size > receiver_size)
    {
        receiver_size=current_receiver_size;
        // OUT(ofs_running, "receiver_size:", receiver_size);
        delete[] receiver_buffer;
        receiver_buffer=new double[receiver_size];
    }

    return 0;
}

void Gint_Gamma::cal_vlocal(
	const double* vlocal_in)
{
	TITLE("Gint_Gamma","cal_vlocal");
	timer::tick("Gint_Gamma","cal_vlocal",'J');

	this->job=cal_local;
	this->vlocal=vlocal_in;
	this->save_atoms_on_grid(GridT);
	if(BFIELD)
	{
		this->gamma_vlocal_B();
	}
	else
	{
		this->gamma_vlocal();
	}

	timer::tick("Gint_Gamma","cal_vlocal",'J');
	return;
}

#include "bfield.h"
void Gint_Gamma::gamma_vlocal_B(void)
{
	TITLE("Grid_Integral","gamma_vlocal_B");

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


	complex<double>** GridVlocal;
	bool perform_gint=true;

	const int lgd_now=GridT.lgd;
	if(lgd_now > 0)
	{
		GridVlocal=new complex<double>*[lgd_now];
		for (int i=0; i<lgd_now; i++)
		{
			GridVlocal[i]=new complex<double>[lgd_now];
			ZEROS(GridVlocal[i], lgd_now);
		}
		Memory::record("Gint_Gamma","GridVlocal",2*lgd_now*lgd_now,"double");
	}
	else if(lgd_now <= 0)
	{
		perform_gint=false;
	}
	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	const double delta_r=ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;

	// allocate 1
	int nnnmax=0;
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax=max(nnnmax, nnn[T]);
	}

	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;
	bool** cal_flag;
	if(max_size!=0)
	{
		dr=new double**[pw.bxyz];
		distance=new double*[pw.bxyz];
		psir_ylm=new double**[pw.bxyz];
		cal_flag=new bool*[pw.bxyz];

		for(int i=0; i<pw.bxyz; i++)
		{
			dr[i]=new double*[max_size];
			distance[i]=new double[max_size];
			psir_ylm[i]=new double*[max_size];
			cal_flag[i]=new bool[max_size];

			ZEROS(distance[i], max_size);
			ZEROS(cal_flag[i], max_size);

			for(int j=0; j<max_size; j++)
			{
				dr[i][j]=new double[3];
				psir_ylm[i][j]=new double[ucell.nwmax];
				ZEROS(dr[i][j],3);
				ZEROS(psir_ylm[i][j],ucell.nwmax);
			}
		}
	}
	double* ylma=new double[nnnmax]; // Ylm for each atom: [bxyz, nnnmax]
	ZEROS(ylma, nnnmax);

	double mt[3]={0,0,0};
	double *vldr3=new double[pw.bxyz];
	double v1=0.0;
	int* vindex=new int[pw.bxyz];
	ZEROS(vldr3, pw.bxyz);
	ZEROS(vindex, pw.bxyz);
	double phi=0.0;

	//This is to store the position of grids, zhiyuan add 2012-01-12//
	double **Rbxyz=new double *[pw.bxyz]; // bxyz*3//
	for(int i=0;i<pw.bxyz;i++)
	{
		Rbxyz[i]=new double[3];
		ZEROS(Rbxyz[i], 3);
	}

	//This is for the phase, zhiyuan add 2012-01-12//
	double A2_ms_A1[3]={0,0,0}; //used to store A2-A1//


	const int nbx=GridT.nbx;
	const int nby=GridT.nby;
	const int nbz_start=GridT.nbzp_start;
	const int nbz=GridT.nbzp;


	for (int i=0; i< nbx; i++)
	{
		if( !perform_gint ) continue;
		for (int j=0; j<nby; j++)
		{
			for (int k=nbz_start; k<nbz_start+nbz; k++) // FFT grid
			{
				this->grid_index=(k-nbz_start) + j * nbz + i * nby * nbz;

				// get the value: how many atoms has orbital value on this grid.
				const int size=GridT.how_many_atoms[ this->grid_index ];
				if(size==0) continue;

				// (1) initialized the phi * Ylm.
				for (int id=0; id<size; id++)
				{
					// there are two parameters we want to know here:
					// in which bigcell of the meshball the atom in?
					// what's the cartesian coordinate of the bigcell?
					const int mcell_index=GridT.bcell_start[grid_index] + id;
					const int imcell=GridT.which_bigcell[mcell_index];

					int iat=GridT.which_atom[mcell_index];

					const int it=ucell.iat2it[ iat ];
					const int ia=ucell.iat2ia[ iat ];

					// meshball_positions should be the bigcell position in meshball
					// to the center of meshball.
					// calculated in cartesian coordinates
					// the vector from the grid which is now being operated to the atom position.
					// in meshball language, is the vector from imcell to the center cel, plus
					// tau_in_bigcell.
					mt[0]=GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
					mt[1]=GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
					mt[2]=GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

					for(int ib=0; ib<pw.bxyz; ib++)
					{
						// meshcell_pos: z is the fastest
						dr[ib][id][0]=GridT.meshcell_pos[ib][0] + mt[0];
						dr[ib][id][1]=GridT.meshcell_pos[ib][1] + mt[1];
						dr[ib][id][2]=GridT.meshcell_pos[ib][2] + mt[2];

						distance[ib][id]=std::sqrt(dr[ib][id][0]*dr[ib][id][0] + dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);
						if(distance[ib][id] <= ORB.Phi[it].getRcut())
						{
							cal_flag[ib][id]=true;
						}
						else
						{
							cal_flag[ib][id]=false;
							continue;
						}

						//if(distance[id] > GridT.orbital_rmax) continue;
						//  Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
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
						const double position=distance[ib][id] / delta_r;
						/*
							this->iq[id]=static_cast<int>(position);
							this->x0[id]=position - static_cast<double>(iq[id]);
							this->x1[id]=1.0 - x0[id];
							this->x2[id]=2.0 - x0[id];
							this->x3[id]=3.0 - x0[id];
							this->x12[id]=x1[id]*x2[id] / 6.0;
							this->x03[id]=x0[id]*x3[id] / 2.0;
						 */
						int ip;
						double dx, dx2, dx3;
						double c1, c2, c3, c4;

						ip=static_cast<int>(position);
						dx=position - ip;
						dx2=dx * dx;
						dx3=dx2 * dx;

						c3=3.0*dx2-2.0*dx3;
						c1=1.0-c3;
						c2=(dx-2.0*dx2+dx3)*delta_r;
						c4=(dx3-dx2)*delta_r;

						//	  int ip=this->iq[id];
						//	  double A=ip+1.0-position/delta_r;
						//	  double B=1.0-A;
						//	  double coef1=(A*A*A-A)/6.0*delta_r*delta_r;
						//	  double coef2=(B*B*B-B)/6.0*delta_r*delta_r;
						Atom* atom1=&ucell.atoms[it];
						for (int iw=0; iw< atom1->nw; iw++)
						{
							if ( atom1->iw2_new[iw] )
							{
								pointer=&ORB.Phi[it].PhiLN(
										atom1->iw2l[iw],
										atom1->iw2n[iw]);
								phi=c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
									+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
							}
							psir_ylm[ib][id][iw]=phi * ylma[atom1->iw2_ylm[iw]];
							//psir_ylm[ib][id][iw]=1;//for test
						}
					}// end ib
				}// end id


				int bindex=0;
				// z is the fastest,
				for(int ii=0; ii<pw.bx; ii++)
				{
					const int iii=i*pw.bx + ii;
					const int ipart=iii*pw.ncy*pw.nczp;
					for(int jj=0; jj<pw.by; jj++)
					{
						const int jjj=j*pw.by + jj;
						const int jpart=jjj*pw.nczp;
						for(int kk=0; kk<pw.bz; kk++)
						{
							const int kkk=k*pw.bz + kk;
							vindex[bindex]=(kkk-pw.nczp_start) + jpart + ipart;
							//						  assert(vindex[bindex] < pw.nrxx);
							++bindex;
						}
					}
				}

				// extract the local potentials.
				for(int ib=0; ib<pw.bxyz; ib++)
				{
					vldr3[ib]=this->vlocal[vindex[ib]] * this->vfactor;
					// vldr3[ib]=1.0e-5; // for test
					// vldr3[bindex]=this->vfactor; // for checking overlap S
				}

				//get the information about position r of grid//
				int index=0;
				double rA[3]={0,0,0};
				for(int ii=0; ii<pw.bx; ii++)
				{
					const double iii=(i*pw.bx + ii)/(double)pw.ncx;
					for(int jj=0; jj<pw.by; jj++)
					{
						const double jjj=(j*pw.by + jj)/(double)pw.ncy;
						for(int kk=0; kk<pw.bz; kk++)
						{
							const double kkk=(k*pw.bz + kk)/(double)pw.ncz;
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
					const int mcell_index1=GridT.bcell_start[grid_index] + ia1;
					const int iat1=GridT.which_atom[mcell_index1];
					const int T1=ucell.iat2it[iat1];
					Atom *atom1=&ucell.atoms[T1];
					const int I1=ucell.iat2ia[iat1];
					// get the start index of local orbitals.
					const int start1=ucell.itiaiw2iwt(T1, I1, 0);

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
						const int mcell_index2=GridT.bcell_start[grid_index] + ia2;
						const int iat2=GridT.which_atom[mcell_index2];
						const int T2=ucell.iat2it[iat2];

						//get phase//
						A2_ms_A1[0]=bfid.A_of_Atom[iat2][0]-bfid.A_of_Atom[iat1][0];
						A2_ms_A1[1]=bfid.A_of_Atom[iat2][1]-bfid.A_of_Atom[iat1][1];
						A2_ms_A1[2]=bfid.A_of_Atom[iat2][2]-bfid.A_of_Atom[iat1][2];
						// only do half part of matrix(including diago part)
						// for T2 > T1, we done all elements, for T2 == T1,
						// we done about half.
						if (T2 >= T1)
						{
							Atom *atom2=&ucell.atoms[T2];
							const int I2=ucell.iat2ia[ GridT.which_atom[mcell_index2]];
							const int start2=ucell.itiaiw2iwt(T2, I2, 0);

							for (int ib=0; ib<pw.bxyz; ib++)
							{
								if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
								{
									double* psi1=psir_ylm[ib][ia1];
									double* psi2=psir_ylm[ib][ia2];
									int iw1_lo=GridT.trace_lo[start1];

									complex<double> phase=0;  //the phase for overlap//
									phase=bfid.cal_phase(A2_ms_A1,Rbxyz[ib],-1);

									// how many orbitals in this type: SZ or DZP or TZP...
									for (int iw=0; iw< atom1->nw; iw++, ++iw1_lo)
									{
										v1=psi1[iw] * vldr3[ib];
										int iw2_lo=GridT.trace_lo[start2];
										complex<double>* result=&GridVlocal[iw1_lo][iw2_lo];
										double* psi2p=psi2;
										double* psi2p_end=psi2p + atom2->nw;
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
											result[0] += phase*v1*psi2p[0];

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
			delete[] cal_flag[i];
		}
		delete[] dr;
		delete[] distance;
		delete[] psir_ylm;
		delete[] cal_flag;
	}

	//delete Rbxyz//
	for(int i=0;i<pw.bxyz;i++)
	{
		delete [] Rbxyz[i];
	}
	delete [] Rbxyz;
	//Parallel to all processors, and upload to H matrix//
	complex<double>* tmp=new complex<double>[NLOCAL];
	for (int i=0; i<NLOCAL; i++)
	{
		ZEROS(tmp, NLOCAL);
		const int mu=GridT.trace_lo[i];
		if (mu >= 0)
		{
			for (int j=0; j<NLOCAL; j++)
			{
				const int nu=GridT.trace_lo[j];
				if (nu >=0)
				{
					if (mu <= nu)
					{
						tmp[j]=GridVlocal[mu][nu];
					}
					else
					{
						//-------------------------------
						// origin:
						// tmp[i]=GridVlocal[nu][mu];
						// mohan fix bug
						// 2011-01-13
						//-------------------------------
						tmp[j] =conj(GridVlocal[nu][mu]);
						//zhiyuan changed 2012-01-13, the v local matrix is Hermite//
					}
				}
			}
		}
		Parallel_Reduce::reduce_complex_double_pool( tmp, NLOCAL);
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}

			// mohan update 2011-04-15
			if(BFIELD)
			{
				LM.set_HSk(i,j,tmp[j],'L');
			}
			else
			{
				//  LM.set_HSgamma(i,j,tmp[j].real(),'L');
			}
		}
	}
	delete[] tmp;

	// mohan update 2010-09-07
	if(GridT.lgd>0)
	{
		for(int i=0; i<GridT.lgd; i++)
		{
			delete[] GridVlocal[i];
		}
		delete[] GridVlocal;
	}


	return;
}




// this subroutine lies in the heart of LCAO algorithms.
// so it should be done very efficiently, very carefully.
// I might repeat again to emphasize this: need to optimize
// this code very efficiently, very carefully.
void Gint_Gamma::gamma_vlocal(void)
{
	TITLE("Gint_Gamma","gamma_vlocal");
	timer::tick("Gint_Gamma","gamma_vlocal",'K');

	bool perform_gint=true;
	double *GridVlocal_pool;
	double **GridVlocal;
	
	//OUT(ofs_running, "start calculate gamma_vlocal");

	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	const double delta_r=ORB.dr_uniform;

	// allocate 1
	int nnnmax=0;
	for(int T=0; T<ucell.ntype; T++)
	{
		nnnmax=max(nnnmax, nnn[T]);
	}

	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	int LD_pool;
	int nblock;
	int *block_size; //band size: number of columns of a band
	int *block_index;
	double *psir_ylm_pool, **psir_ylm;
	double *psir_vlbr3_pool, **psir_vlbr3;
    int *at;
    int **cal_flag;
	
	double* ylma;

	double *vldr3;
	double v1=0.0;
	int* vindex;

	const int nbx=GridT.nbx;
	const int nby=GridT.nby;
	const int nbz_start=GridT.nbzp_start;
	const int nbz=GridT.nbzp;

	const int ncyz=pw.ncy*pw.nczp;

	const int lgd_now=GridT.lgd;	
	if(max_size<=0 || lgd_now <= 0) 
	{
		perform_gint=false;
		goto ENDandRETURN;
	}
	GridVlocal_pool=new double [lgd_now*lgd_now];
	ZEROS(GridVlocal_pool, lgd_now*lgd_now);
	GridVlocal=new double*[lgd_now];
	for (int i=0; i<lgd_now; i++)
	{
		GridVlocal[i]=&GridVlocal_pool[i*lgd_now];
	}
	Memory::record("Gint_Gamma","GridVlocal",lgd_now*lgd_now,"double");
	
	 ylma=new double[nnnmax]; // Ylm for each atom: [bxyz, nnnmax]
	ZEROS(ylma, nnnmax);
	vldr3=new double[pw.bxyz];
	vindex=new int[pw.bxyz];
	ZEROS(vldr3, pw.bxyz);
	ZEROS(vindex, pw.bxyz);
	
	LD_pool=max_size*ucell.nwmax;
	distance=new double*[pw.bxyz];
    block_size=new int[max_size];
    block_index=new int[max_size+1];
    at=new int[max_size];
	psir_ylm_pool=new double[pw.bxyz*LD_pool];
	psir_ylm=new double *[pw.bxyz];
	ZEROS(psir_ylm_pool, pw.bxyz*LD_pool);
	psir_vlbr3_pool=new double[pw.bxyz*LD_pool];
	psir_vlbr3=new double *[pw.bxyz];
	ZEROS(psir_vlbr3_pool, pw.bxyz*LD_pool);
    cal_flag=new int*[pw.bxyz];
	for(int i=0; i<pw.bxyz; i++)
	{
		distance[i]=new double[max_size];
		psir_ylm[i]=&psir_ylm_pool[i*LD_pool];
		psir_vlbr3[i]=&psir_vlbr3_pool[i*LD_pool];
		ZEROS(distance[i], max_size);
        cal_flag[i]=new int[max_size];
	}

	int *block_iw; // index of wave functions of each block;
	block_iw=new int[max_size];
		
	for (int i=0; i< nbx; i++)
	{
		const int ibx=i*pw.bx;
		for (int j=0; j<nby; j++)
		{
			const int jby=j*pw.by; 
			for (int k=nbz_start; k<nbz_start+nbz; k++) // FFT grid
			{
				//OUT(ofs_running, "====================");
				//OUT(ofs_running, "i", i);
				//OUT(ofs_running, "j", j);
				//OUT(ofs_running, "k", k);
				this->grid_index=(k-nbz_start) + j * nbz + i * nby * nbz;

				// get the value: how many atoms has orbital value on this grid.
				const int size=GridT.how_many_atoms[ this->grid_index ];
				if(size==0) continue;
				const int kbz=k*pw.bz-pw.nczp_start;
				setVindex(ncyz, ibx, jby, kbz, vindex);
				//OUT(ofs_running, "vindex was set");
				// extract the local potentials.
				for(int ib=0; ib<pw.bxyz; ib++)
				{
					vldr3[ib]=this->vlocal[vindex[ib]] * this->vfactor;
				}
				
				//OUT(ofs_running, "vldr3 was inited");
                timer::tick("Gint_Gamma","psir_vlocal",'J');
                cal_psir_ylm(size, this->grid_index, delta_r, distance, ylma,
                        at, block_index, block_iw, block_size, 
                        cal_flag, psir_ylm);
                timer::tick("Gint_Gamma","psir_vlocal",'J');
				timer::tick("Gint_Gamma","meshball_vlocal",'J');
				cal_meshball_vlocal(size, LD_pool, block_iw, block_size, block_index, cal_flag, vldr3, psir_ylm, psir_vlbr3, vindex, lgd_now, GridVlocal);
				timer::tick("Gint_Gamma","meshball_vlocal",'J');
				//OUT(ofs_running, "GridVlocal was calculated");
			}// k
		}// j
	}// i
	
	for(int i=0; i<pw.bxyz; i++)
	{
        delete[] cal_flag[i];
		delete[] distance[i];
	}
    delete[] cal_flag;
	delete[] vindex;
	delete[] ylma;
	delete[] vldr3;	
	delete[] block_iw;
	delete[] distance;
	delete[] psir_vlbr3;
	delete[] psir_vlbr3_pool;
	delete[] psir_ylm;
	delete[] psir_ylm_pool;
    delete[] block_index;
    delete[] block_size;
    delete[] at;
	//OUT(ofs_running, "temp variables are deleted");

ENDandRETURN:    
    timer::tick("Gint_Gamma","gamma_vlocal",'K');
    MPI_Barrier(MPI_COMM_WORLD);
    timer::tick("Gint_Gamma","distri_vl",'K');

    // setup send buffer and receive buffer size
    // OUT(ofs_running, "Start transforming vlocal from grid distribute to 2D block");
    if(chr.new_e_iteration)
    {
        timer::tick("Gint_Gamma","distri_vl_index",'K');
        // OUT(ofs_running, "Setup Buffer Parameters");
        // inline int setBufferParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk,
        //                             int& sender_index_size, int*& sender_local_index, 
        //                             int*& sender_size_process, int*& sender_displacement_process, 
        //                             int& sender_size, double*& sender_buffer,
        //                             int& receiver_index_size, int*& receiver_global_index, 
        //                             int*& receiver_size_process, int*& receiver_displacement_process, 
        //                             int& receiver_size, double*& receiver_buffer)
        setBufferParameter(ParaO.comm_2D, ParaO.blacs_ctxt, ParaO.nb,
                           ParaO.sender_index_size, ParaO.sender_local_index, 
                           ParaO.sender_size_process, ParaO.sender_displacement_process, 
                           ParaO.sender_size, ParaO.sender_buffer,
                           ParaO.receiver_index_size, ParaO.receiver_global_index, 
                           ParaO.receiver_size_process, ParaO.receiver_displacement_process, 
                           ParaO.receiver_size, ParaO.receiver_buffer);
        OUT(ofs_running, "vlocal exchange index is built");
        OUT(ofs_running, "buffer size(M):", (ParaO.sender_size+ParaO.receiver_size)*sizeof(double)/1024/1024);
        OUT(ofs_running, "buffer index size(M):", (ParaO.sender_index_size+ParaO.receiver_index_size)*sizeof(int)/1024/1024);
        timer::tick("Gint_Gamma","distri_vl_index",'K');
    }

    // OUT(ofs_running, "Start data transforming");
    timer::tick("Gint_Gamma","distri_vl_value",'K');
    // put data to send buffer
    for(int i=0; i<ParaO.sender_index_size; i+=2)
    {
        // const int idxGrid=ParaO.sender_local_index[i];
        // const int icol=idxGrid%lgd_now; 
        // const int irow=(idxGrid-icol)/lgd_now;
        const int irow=ParaO.sender_local_index[i];
        const int icol=ParaO.sender_local_index[i+1];
        if(irow<=icol)
            ParaO.sender_buffer[i/2]=GridVlocal[irow][icol];
        else
            ParaO.sender_buffer[i/2]=GridVlocal[icol][irow];
    }
    // OUT(ofs_running, "vlocal data are put in sender_buffer, size(M):", ParaO.sender_size*8/1024/1024);

    // use mpi_alltoall to get local data
    MPI_Alltoallv(ParaO.sender_buffer, ParaO.sender_size_process, ParaO.sender_displacement_process, MPI_DOUBLE, 
                  ParaO.receiver_buffer, ParaO.receiver_size_process, ParaO.receiver_displacement_process, MPI_DOUBLE, ParaO.comm_2D);
    // OUT(ofs_running, "vlocal data are exchanged, received size(M):", ParaO.receiver_size*8/1024/1024);
    // put local data to H matrix
    for(int i=0; i<ParaO.receiver_index_size; i+=2)
    {
        // int g_col=ParaO.receiver_global_index[i]%NLOCAL;
        // int g_row=(ParaO.receiver_global_index[i]-g_col)/NLOCAL;
        const int g_row=ParaO.receiver_global_index[i];
        const int g_col=ParaO.receiver_global_index[i+1];
        // if(g_col<0 || g_col>=NLOCAL||g_row<0 || g_row>=NLOCAL) 
        // {
        //     OUT(ofs_running, "index error, i:", i);
        //     OUT(ofs_running, "index：", ParaO.receiver_global_index[i]);
        //     OUT(ofs_running, "g_col:", g_col);
        //     OUT(ofs_running, "g_col:", g_col);
        // }
        LM.set_HSgamma(g_row,g_col,ParaO.receiver_buffer[i/2],'L');
    }
    // OUT(ofs_running, "received vlocal data are put in to H")
    timer::tick("Gint_Gamma","distri_vl_value",'K');
    timer::tick("Gint_Gamma","distri_vl",'K');

    //OUT(ofs_running, "reduce all vlocal ok,");
    if(perform_gint)
    {
        delete[] GridVlocal_pool;
        delete[] GridVlocal;
    }
    //OUT(ofs_running, "ALL GridVlocal was calculated");
    return;
}
