#include "subgrid_oper.h"
#include "../src_pw/tools.h"
#include "../src_pw/global.h"

SubGrid_oper::SubGrid_oper()
{
	trace_lo_tot = new int[1];	
	lgd=0;
	allocate_totwfc = false;
}
SubGrid_oper::~SubGrid_oper()
{
	delete[] trace_lo_tot;
}

//--------------------------------------------------
// because only the DIAG_WORLD processors
// have augmented wave functions,
// others will not have, we use this 
// augmented wave functions to calculate
// the force, not related to the grid integration
// Only DIAG_WORLD use density matrix (2D) only.
//--------------------------------------------------
/*
void SubGrid_oper::cal_totwfc_aug()
{
	TITLE("SubGrid_oper","cal_totwfc_aug");

	int *occupy = new int[NLOCAL];
	ZEROS(occupy, NLOCAL);
	for(int i=0; i<NLOCAL; ++i)
	{
		if(LOWF.trace_aug[i]>=0)
		{
			occupy[i] += 1;
		}
	}

	// reduce occupy and get the full occupations.
#ifdef __MPI
	Parallel_Reduce::reduce_int_grid(occupy, NLOCAL);
#endif

	delete[] trace_aug_tot;
	trace_aug_tot = new int[NLOCAL];
	ZEROS(trace_aug_tot,NLOCAL);
	for(int i=0; i<NLOCAL; ++i)
	{
		trace_aug_tot[i]=LOWF.trace_aug[i];
	}

	// how many occupied wave functions.	
#ifdef __MPI
	Parallel_Reduce::reduce_int_grid(trace_aug_tot, NLOCAL);
#endif

	if(GRANK==0)
	{
		// rescale trace_lo_tot
		for(int i=0; i<NLOCAL; ++i)
		{
			if(occupy[i]!=0)
			{
				trace_aug_tot[i]/=occupy[i];
			}
		}	
		delete[] occupy;

		ofs_running << " Second" << endl;
		for(int i=0; i<NLOCAL; ++i)
		{
			ofs_running << " i=" << i << " trace_lo_tot=" << trace_lo_tot[i] << endl;
		}

		//-----------------------------------------
		// calculate the total lgd.
		//-----------------------------------------
		this->lgd_aug = 0;
		for(int i=0; i<NLOCAL; ++i)
		{
			if(trace_aug_tot[i]>=0)
			{
				++lgd_aug;
			}
		}

		OUT(ofs_running,"SGO.lgd_aug",lgd);
	


	return;
}
*/

void SubGrid_oper::cal_totwfc()
{
	TITLE("SubGrid_oper","cal_totwfc");

	//-----------------------------------------
	// combine the wave functions index in the 
	// 'Grid group'
	//-----------------------------------------
	int *occupy = new int[NLOCAL];
	ZEROS(occupy, NLOCAL);	
	for(int iw=0; iw<NLOCAL; ++iw)
	{
		if(GridT.trace_lo[iw] >= 0)
		{
			occupy[iw] = 1;
		}
	}		

	// reduce occupy and get the full occupations.
#ifdef __MPI
	// mohan 2012-02-23
	Parallel_Reduce::reduce_int_grid(occupy, NLOCAL);
#endif

	/*
	for(int i=0; i<NLOCAL; ++i)
	{
		ofs_running << " i=" << i << " occupy=" << occupy[i] << endl;
	}
	*/

	delete[] trace_lo_tot;
	trace_lo_tot = new int[NLOCAL];
	ZEROS(trace_lo_tot, NLOCAL);


	int count = 0;
	for(int iw=0; iw<NLOCAL; ++iw)
	{
		// mohan fix bug 2012-02-23, shoud > 0,
		// not >= 0.
		if(occupy[iw]>0)
		{
			trace_lo_tot[iw]=count;
			++count;
		}
		else
		{
			trace_lo_tot[iw]=-1;
		}
	}
	//-----------------------------------------
	// calculate the total lgd.
	//-----------------------------------------
	this->lgd = count;
	
	//------------
	// for test
	//------------
	/*
	ofs_running << " trace_lo_tot" << endl;
	for(int iw=0; iw<NLOCAL; ++iw)
	{
		ofs_running << " iw=" << iw << " trace_lo_tot=" << trace_lo_tot[iw] 
		<< " trace_lo=" << GridT.trace_lo[iw] 
		<< " occupy=" << occupy[iw] 
		<< endl;
	}
	*/

	if(GRANK==0)
	{
          	//xiaohui add 'OUT_LEVEL', 2015-09-16
		if(OUT_LEVEL != "m") OUT(ofs_running,"SGO.lgd",lgd);
	
		if(lgd==0)
		{
			// need to allocate the pointer,
			// because it would be some array's
			// entrance parameters.
			// mohan update 2021-02-12
			this->totwfc = new double**[1];
			for(int is=0; is<1; ++is)
			{
				this->totwfc[is] = new double*[NBANDS];
			}
			this->allocate_totwfc=false;
			return;
		}
	
		assert(this->allocate_totwfc==false);

		static bool first_time=true;
		if(first_time)
		{
//			LOWF.init_Cij(1);
			first_time=false;
		}


		this->totwfc = new double**[1];
		for(int is=0; is<1; ++is)
		{
			this->totwfc[is] = new double*[NBANDS];
			for(int ib=0; ib<NBANDS; ++ib)
			{
				this->totwfc[is][ib] = new double[lgd];	
				ZEROS(totwfc[is][ib], lgd);

				// mohan update 2012-02-10
				//if(DIAGO_TYPE!="cg") xiaohui modify 2013-09-02
				if(KS_SOLVER!="cg") //xiaohui add 2013-09-02
				{	
					for(int i=0; i<NLOCAL; ++i)
					{
						if(occupy[i])
						{
							if(ib==i)
							{
								totwfc[is][ib][trace_lo_tot[i]] = 1.0;
							}
						}
					}
				}
				else //use cg method
				{
					// bug: dimension of WFC_GAMMA is smaller than totwfc
					for(int i=0; i<lgd; ++i)
					{
		// mohan comment out 2021-02-09
		//				this->totwfc[is][ib][i] = LOWF.WFC_GAMMA[CURRENT_SPIN][ib][i]; //mohan update 2012-02-07	
					}
				}
			}	
		}

		allocate_totwfc = true;
	}

	delete[] occupy; //mohan fix bug 2012-03-25

	return;
}


void SubGrid_oper::dis_subwfc()
{
	TITLE("SubGrid_oper","dis_subwfc");

#ifdef __MPI

//	cout << " distribute the wave functions " << endl;

	//------------------------------------------
	// bcast the eigenvalues
	//------------------------------------------
	for(int ik=0; ik<kv.nks; ++ik)
	{
		MPI_Bcast(wf.ekb[ik], NBANDS, MPI_DOUBLE, 0, GRID_WORLD);
	}	

	MPI_Status status;

	for(int i=0; i<GSIZE; ++i)
	{
		if(GRANK==0)
		{
			if(i==0)
			{
				//---------------------------------------------
				// Transfer the data from totwfc to WFC_GAMMA.
				//---------------------------------------------
				for(int iw=0; iw<NLOCAL; ++iw)
				{
					const int mu1 = GridT.trace_lo[iw];
					if(mu1 >= 0)
					{
						const int mu2 = this->trace_lo_tot[iw];

						for(int ib=0; ib<NBANDS; ++ib)
						{
		// mohan comment out 2021-02-09
		//					LOWF.WFC_GAMMA[CURRENT_SPIN][ib][mu1] = this->totwfc[0][ib][mu2];
						}//ib
					}//mu1>=0
				}//iw
			}//i
			else
			{
				int tag;
				// receive trace_lo2
				tag = i * 10;
				int* trace_lo2 = new int[NLOCAL];
				MPI_Recv(trace_lo2, NLOCAL, MPI_INT, i, tag, GRID_WORLD, &status);

/*
				ofs_running << " Proc " << i << endl;
				for(int i=0; i<NLOCAL; ++i)
				{
					ofs_running << setw(5) << i << setw(10) << trace_lo2[i] << endl;
				}
				*/

				// receive lgd2
				int lgd2 = 0;
				tag = i * 10 + 1;
				MPI_Recv(&lgd2, 1, MPI_INT, i, tag, GRID_WORLD, &status);

//				ofs_running << " receive=" << tag << endl;
				
				// send csend
				double* csend = new double[NBANDS*lgd2];
				ZEROS(csend, NBANDS*lgd2);

				for (int iw=0; iw<NLOCAL; iw++)
				{
					const int mu1 = trace_lo2[iw];
					if (mu1>=0)
					{
						const int mu2 = this->trace_lo_tot[iw]; 
						if(mu2<0)
						{
							ofs_running << " iw = " << iw << endl;
							ofs_running << " GridT.trace_lo=" << mu1 << endl;
							ofs_running << " trace_lo_tot=" << mu2 << endl;
							assert(mu2>=0);
						}

						for (int ib=0; ib<NBANDS; ib++)
						{
							csend[mu1*NBANDS+ib] = this->totwfc[0][ib][mu2];
						}
					}
				}

				tag = i * 10 + 2;
//				ofs_running << " send=" << tag << endl;
				MPI_Send(csend,NBANDS*lgd2,MPI_DOUBLE,i,tag,GRID_WORLD);
//				ofs_running << " send done." << endl;

				delete[] csend;
				delete[] trace_lo2;
			}
		}//GRANK=0
		else if(i==GRANK)
		{
			int tag;
			// send trace_lo
			tag = GRANK * 10;
			MPI_Send(GridT.trace_lo, NLOCAL, MPI_INT, 0, tag, GRID_WORLD);

			//ofs_running << " send1." << endl;

			// send GridT.lgd
			tag = GRANK * 10 + 1;
			MPI_Send(&GridT.lgd, 1, MPI_INT, 0, tag, GRID_WORLD);

			//ofs_running << " send2." << endl;

			// receive c
			double* crecv = new double[NBANDS*GridT.lgd];
			ZEROS(crecv, NBANDS*GridT.lgd);

			tag = GRANK * 10 + 2;
//			ofs_running << " receive=" << tag << endl;
			MPI_Recv(crecv, NBANDS*GridT.lgd, MPI_DOUBLE, 0, tag, GRID_WORLD, &status);
//			ofs_running << " receive done." << endl;


			for(int ib=0; ib<NBANDS; ++ib)
			{
				for (int mu=0; mu<GridT.lgd; ++mu)
				{
		// mohan comment out 2021-02-09
		//			LOWF.WFC_GAMMA[CURRENT_SPIN][ib][mu] = crecv[mu*NBANDS+ib];
				}
			}

			delete[] crecv;
		}
		MPI_Barrier(GRID_WORLD);
	}//end i


	//-------------------
	// Test
	//-------------------
	/*
	ofs_running << " WFC " << " CURRENT_SPIN=" << CURRENT_SPIN << endl;
	for(int i=0; i<NBANDS; ++i)
	{
		for(int j=0; j<NLOCAL; ++j)
		{
			const int mu = GridT.trace_lo[j];
			if(mu>=0)
			{
				//				if( abs(LOWF.WFC_GAMMA[0][i][mu] > 1.0e-8) )
				ofs_running << setw(5) << i+1 << setw(8) << j+1 
					<< setw(15) << LOWF.WFC_GAMMA[CURRENT_SPIN][i][mu] << endl; 
			}
		}
	}
	*/

#else

// this is for serial version

	//---------------------------------------------
	// Transfer the data from totwfc to WFC_GAMMA.
	//---------------------------------------------
	for(int iw=0; iw<NLOCAL; ++iw)
	{
		const int mu1 = GridT.trace_lo[iw];
		if(mu1 >= 0)
		{
			const int mu2 = this->trace_lo_tot[iw];

			for(int ib=0; ib<NBANDS; ++ib)
			{
		// mohan comment out 2021-02-09
//				LOWF.WFC_GAMMA[CURRENT_SPIN][ib][mu1] = this->totwfc[0][ib][mu2];
			}//ib
		}//mu1>=0
	}//iw


#endif //mohan fix bug 2012-02-04

	if(GRANK==0)
	{
		// 1 spin because it can not be calculated two spin
		// at the same time.
		for(int is=0; is<1; ++is)
		{
			for(int ib=0; ib<NBANDS; ++ib)
			{
				if(allocate_totwfc)
				{
					delete[] this->totwfc[is][ib];
				}
			}
			delete[] this->totwfc[is];
		}
		delete[] this->totwfc;
		this->allocate_totwfc = false;
	}
	return;
}

