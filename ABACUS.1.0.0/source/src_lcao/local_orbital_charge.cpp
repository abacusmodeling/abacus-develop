#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "blas_interface.h"

//#include "../src_onscaling/on_tests.h"
//#include "../src_siao/selinv.h"
// 2014.10.29 add memory pool for DM and DM_B by yshen
Local_Orbital_Charge::Local_Orbital_Charge()
{
	// for gamma algorithms.
	this->init_DM = false;	
	this->lgd_now = 0;
	this->lgd_last = 0;

	// for k-dependent algorithms.
	this->init_DM_R = false;
	out_dm = 0;
	//xiaohui add 2014-06-19
	//local_band = new int[1];
	//Z_wg = new double[1];
	//Z_LOC = new double[1];
}

Local_Orbital_Charge::~Local_Orbital_Charge()
{
	// with gamma point only
	 if (this->init_DM)
	 {
		if(BFIELD)
		{
			for (int is=0; is<NSPIN; is++)
			{
				delete[] DM_B[is];
				delete[] DM_B_pool[is];
			}
			delete[] DM_B;
			delete[] DM_B_pool;
		}
		else
		{
			for (int is=0; is<NSPIN; is++)
			{
				delete[] DM[is];
				delete[] DM_pool[is];
			}
			delete[] DM;
			delete[] DM_pool;
		}
	 }

	// with k points
	if (this->init_DM_R)
	{
		for(int is=0; is<NSPIN; is++)
		{
			delete[] DM_R[is];
		}
		delete[] DM_R;
	}

	//xiaohui add 2014-06-19
	//delete[] loc_bands;
	//delete[] Z_wg;
	//delete[] Z_LOC;
/*
	//xiaohui add 2014-06-20
	for (int is=0; is<NSPIN; is++)
	{
		for (int i=0; i<NLOCAL; i++)
		{
			delete[] DM[is][i]; // mohan fix 2009-08-20
		}
		delete[] DM[is];
	}
	delete[] DM;
*/
}

#include "lcao_nnr.h"
void Local_Orbital_Charge::allocate_DM_k(void)
{
	 TITLE("Local_Orbital_Charge","allocate_k");

	this->nnrg_now = LNNR.nnrg;
	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"nnrg_last",nnrg_last);
	if(OUT_LEVEL != "m") OUT(ofs_running,"nnrg_now",nnrg_now);

	if(this->init_DM_R)
	{
		assert(nnrg_last > 0);
		for(int is=0; is<NSPIN; is++)
		{
			delete[] DM_R[is];
		}
		delete[] DM_R;
		init_DM_R=false;
	}

	if(nnrg_now>0)
	{
		this->DM_R = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			this->DM_R[is] = new double[nnrg_now];
			ZEROS(DM_R[is], nnrg_now);
		}
		this->nnrg_last = nnrg_now;
		this->init_DM_R = true;
		Memory::record("LocalOrbital_Charge","Density_Matrix",NSPIN*nnrg_now,"double");
	}
	else if(nnrg_now==0)
	{
		this->init_DM_R = false;
	}
	else
	{
		WARNING_QUIT("Local_Orbital_Charge::allocate_k","check init_DM_R.");
	}

	return;
}

// allocate density kernel may change once the ion
// positions change
void Local_Orbital_Charge::allocate_gamma(const Grid_Technique &gt)
{
	 TITLE("Local_Orbital_Charge","allocate_gamma");

	//xiaohui add 2014-06-20
	//for(int is=0; is<NSPIN; is++)
	//{
	//if(!this->init_DM)
	//{
	//	this->DM = new double**[NSPIN];
	//	for(int is=0; is<NSPIN; is++)
	//	{
	//		this->DM[is] = new double*[NLOCAL];

	//		for (int i=0; i<NLOCAL; i++)
	//		{
	//			DM[is][i] = new double[NLOCAL];
	//			ZEROS(DM[is][i], NLOCAL);
	//		}

	//		this->init_DM = true;
	//		Memory::record("LocalOrbital_Charge","Density_Kernal",NSPIN*NLOCAL*NLOCAL,"double");
	//	}
	//}

		//for (int i=0; i<NLOCAL; i++)
		//{
		//	ZEROS(this->DM[is][NLOCAL], NLOCAL);
		//}
	//}

	//xiaohui modify 2014-06-18

	// mohan fix serious bug 2010-09-06
	this->lgd_now = gt.lgd;
	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"lgd_last",lgd_last);
	if(OUT_LEVEL != "m") OUT(ofs_running,"lgd_now",lgd_now);

	// mohan add 2010-07-01
	if(this->init_DM)
	{
		if(BFIELD)
		{
			assert(lgd_last > 0);
			for (int is=0; is<NSPIN; is++)
			{
				delete[] DM_B[is];
			    delete[] DM_B_pool[is];
			}
			delete[] DM_B;
			delete[] DM_B_pool;
		}
		else
		{
			assert(lgd_last > 0);
			for (int is=0; is<NSPIN; is++)
			{
				delete[] DM[is];
			    delete[] DM_pool[is];
			}
			delete[] DM;
			delete[] DM_pool;
		}
		init_DM = false;
	}

	 assert(lgd_now <= NLOCAL);

	// mohan update 2010-09-06
	if(lgd_now > 0)
	{
		if(BFIELD)
		{
			this->DM_B = new complex<double> **[NSPIN];
			this->DM_B_pool = new complex<double> *[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				this->DM_B_pool[is]=new complex<double> [lgd_now*lgd_now];
				ZEROS(DM_B_pool[is], lgd_now*lgd_now);
				 this->DM_B[is] = new complex<double>*[lgd_now];

				 for (int i=0; i<lgd_now; i++)
				 {
					DM_B[is][i] = &DM_B_pool[is][i*lgd_now];
				 }
				 Memory::record("LocalOrbital_Charge","Density_Kernal",NSPIN*lgd_now*lgd_now,"complex<double>");
			}
		}
		else
		{
			this->DM = new double**[NSPIN];
			this->DM_pool = new double *[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				this->DM_pool[is]=new double [lgd_now*lgd_now];
				ZEROS(DM_pool[is], lgd_now*lgd_now);
				this->DM[is] = new double*[lgd_now];

				for (int i=0; i<lgd_now; i++)
				{
					DM[is][i] = &DM_pool[is][i*lgd_now];
				}
				Memory::record("LocalOrbital_Charge","Density_Kernal",NSPIN*lgd_now*lgd_now,"double");
			}
		}
		this->init_DM = true;
		this->lgd_last = lgd_now;
		//xiaohui add 'OUT_LEVEL', 2015-09-16
		if(OUT_LEVEL != "m") ofs_running << " allocate DM , the dimension is " << lgd_now << endl;
	}
	else if(lgd_now == 0)
	{
		this->init_DM = false;
	}
	else
	{
		WARNING_QUIT("Local_Orbital_Charge::allocate","lgd<0!Something Wrong!");
	}
	 return;
}

void Local_Orbital_Charge::sum_bands(void)
{
	TITLE("Local_Orbital_Charge","sum_bands");
	timer::tick("Local_Orbital_Cha","sum_bands",'E');

	en.eband = 0.0;
	//xiaohui modify 2013-09-02
	//if(LINEAR_SCALING == 2)
	//{
	//	//en.eband = ON.eband;
	//} 
	//else if(LINEAR_SCALING == 1)
	//{
	//	//xiaohui modified 2013-07-22
	//	//if(DIAGO_TYPE=="selinv")
	//	//{
	//	//	en.eband = Selinv::eband;
	//	//}
	//	//else
	//	//{
	//	//	ofs_running << " calculate eband " << endl;
	//		for(int ik=0; ik<kv.nks; ik++)
	//		{
	//			for (int ib=0; ib<NBANDS; ib++)
	//			{
	//				en.eband += wf.ekb[ik][ib] * wf.wg(ik, ib);
	//			//ofs_running << setw(15) << wf.ekb[ik][ib] << setw(15) << wf.wg(ik,ib) << endl; 
	//			}//ib
	//		}//ik
	//	//}
	//}
	//else
	//{
	//	WARNING_QUIT("Local_Orbital_Charge","check the parameter: LINEAR SCALINIG");
	//} xiaohui modify 2013-09-02. Attention...

	//xiaohui add 2013-09-02
	for(int ik=0; ik<kv.nks; ik++)
	{
		for (int ib=0; ib<NBANDS; ib++)
		{
			en.eband += wf.ekb[ik][ib] * wf.wg(ik, ib);
		}
	} //xiaohui add 2013-09-02. Attention...


	//------------------------------------------------------------
	 //calculate density matrix, using coefficients of local basis.
	//------------------------------------------------------------
	//xiaohui modify 2013-09-02
	//if(LINEAR_SCALING == 2)
	//{
	//	// density matrix has already been calculated.
	//}
	//else
	//{
	//	if(GAMMA_ONLY_LOCAL)
	//	{
	//		if(DIAGO_TYPE=="selinv")
	//		{
	//			//density matrix has already been calcualted.
	//		}
	//		else
	//		{
	//			this->cal_dk_gamma();//calculate the density matrix.
	//		}
	//		// @@@@@@@
	//		// test
	//		// @@@@@@@
	//		/*
	//		cout << " Density Matrix:";
	//		double sum = 0.0;
	//		for(int i=0; i<NLOCAL; i++)
	//		{
	//			cout << endl;
	//			for(int j=0; j<NLOCAL; j++)
	//			{
	//				cout << setw(15) << DM[0][i][j];
	//				sum += DM[0][i][j];
	//			}
	//		}
	//		cout << "\n Sum of density kernel = " << sum << endl;
	//		*/
	//	}
	//	else
	//	{
	//		NOTE("Calculate the density matrix!");
	//		this->cal_dk_k( GridT );
	//	}
	//} xiaohui modify 2013-09-02. Attention...

	//xiaohui add 2013-09-02
	if(GAMMA_ONLY_LOCAL)
	{
		if(KS_SOLVER=="selinv")
		{
			//density matrix has already been calcualted.
		}
		else
		{
			//xiaohui modify 2014-06-18
			this->cal_dk_gamma();//calculate the density matrix.
		}
	}
	else
	{
		NOTE("Calculate the density matrix!");
		this->cal_dk_k( GridT );
	} //xiaohui add 2013-09-02. Attention...
	 
	for(int is=0; is<NSPIN; is++)
	{
	 	ZEROS( chr.rho[is], pw.nrxx ); // mohan 2009-11-10
	}

	//------------------------------------------------------------
	//calculate the density in real space grid.
	//------------------------------------------------------------
	 time_t start = time(NULL);

	if(GAMMA_ONLY_LOCAL)
	{
		if(GRID_SPEED==1)
		{
				UHM.GG.cal_rho();
		}
		else if(GRID_SPEED==2)
		{
			UHM.GS.cal_rho();
		}
	}
	else
	{
		NOTE("Calculate the charge on real space grid!");
		UHM.GK.calculate_charge();
	}

	 time_t end = time(NULL);

	 //ofs_running << " START_Charge Time : " << ctime(&time_charge_start);
	 //ofs_running << " END_Charge	Time : " << ctime(&time_charge_end);
	 //ofs_running << " FINAL_Charge Time : " << difftime(time_charge_end, time_charge_start) << " (Seconds)" << endl;

	OUT_TIME("charge grid integration", start, end);

	//BLOCK_HERE("sum_bands::before renormalize rho");	

	 chr.renormalize_rho();

	timer::tick("Local_Orbital_Cha","sum_bands",'E');
	 return;
}

#include "record_adj.h"

inline void cal_DM_ATOM(const Grid_Technique &gt, const complex<double> fac, Record_adj RA,
				   const int ia1, const int iw1_lo, const int nw1, const int gstart, 
				   complex<double> *WFC_PHASE, complex<double> **DM_ATOM)
{
	const char transa='N', transb='T';	
	const complex<double> alpha=1, beta=1;
	int ispin=0;

	for(int ik=0; ik<kv.nks; ik++)
	{
		complex<double> **wfc = LOWF.WFC_K[ik];
		if(NSPIN==2) 
		{
			ispin = kv.isk[ik];
		}
		int atom2start=0;

		for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
		{
			complex<double> *DM=&DM_ATOM[ispin][atom2start];
			const int T2 = RA.info[ia1][ia2][3];
			const int I2 = RA.info[ia1][ia2][4];
			Atom* atom2 = &ucell.atoms[T2];
			const int start2 = ucell.itiaiw2iwt(T2,I2,0);
			const int iw2_lo=gt.trace_lo[start2];
			const int nw2=atom2->nw;
			complex<double> exp_R= exp( fac * (
						kv.kvec_d[ik].x * RA.info[ia1][ia2][0] + 
						kv.kvec_d[ik].y * RA.info[ia1][ia2][1] + 
						kv.kvec_d[ik].z * RA.info[ia1][ia2][2]  
						) );
			
			//ZEROS(WFC_PHASE, NBANDS*nw1);
			int ibStart=0;
			int nRow=0;
			for(int ib=0; ib<NBANDS; ++ib)
			{
				const double w1=wf.wg(ik,ib);
				if(w1>0)
				{
					if(nRow==0) ibStart=ib;
					const int iline=nRow*nw1;
					complex<double> phase=exp_R*w1;
					for(int iw1=0; iw1<nw1; ++iw1)
						WFC_PHASE[iline+iw1]=phase*conj(wfc[ib][iw1_lo+iw1]);
					++nRow;
				}
				else
					break;
			} // ib
			zgemm_(&transa, &transb, &nw2, &nw1, &nRow, &alpha,
				&wfc[ibStart][iw2_lo], &gt.lgd, 
				WFC_PHASE, &nw1,
				&beta, DM, &nw2);			
			atom2start+=nw1*nw2;
		} // ia2
	} // ik
	return;
}


void Local_Orbital_Charge::cal_dk_k(const Grid_Technique &gt)
{
	TITLE("Local_Orbital_Charge","cal_dk_k");
	timer::tick("LCAO_Charge","cal_dk_k",'F');	
	int nnrg = 0;
	Vector3<double> tau1, dtau;
		
	Record_adj RA;
	RA.for_grid(gt);

	int ca = 0;
	complex<double> fac = TWO_PI * IMAG_UNIT;

	complex<double> *WFC_PHASE=new complex<double>[NLOCAL*ucell.nwmax];
	
	int DM_ATOM_SIZE=1;	
	complex<double> **DM_ATOM=new complex<double> *[NSPIN];
	for(int is=0; is<NSPIN; ++is)
	{
		DM_ATOM[is]=new complex<double>[DM_ATOM_SIZE];
		ZEROS(DM_ATOM[is], DM_ATOM_SIZE);
	}
	for(int T1=0; T1<ucell.ntype; T1++)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for(int I1=0; I1<atom1->na; I1++)
		{
			const int iat = ucell.itia2iat(T1,I1);
			if(gt.in_this_processor[iat])
			{
				const int start1 = ucell.itiaiw2iwt(T1,I1,0);
				const int gstart = LNNR.nlocstartg[iat];
				const int ng = LNNR.nlocdimg[iat];
				const int iw1_lo=gt.trace_lo[start1];
				const int nw1=atom1->nw;

				if(DM_ATOM_SIZE<ng)
				{
					DM_ATOM_SIZE=ng;
					for(int is=0; is<NSPIN; ++is)
						delete[] DM_ATOM[is];
					for(int is=0; is<NSPIN; ++is)
						DM_ATOM[is]=new complex<double>[DM_ATOM_SIZE];
				}
				for(int is=0; is<NSPIN; ++is)
					ZEROS(DM_ATOM[is], ng);
				ZEROS(WFC_PHASE, NBANDS*nw1);
				cal_DM_ATOM(gt, fac, RA, ca, iw1_lo, nw1, gstart, WFC_PHASE, DM_ATOM);

				++ca;

				for(int is=0; is<NSPIN; ++is)
				{
					for(int iv=0; iv<ng; ++iv)
					{
						this->DM_R[is][gstart+iv]=DM_ATOM[is][iv].real();
					}
				}
			} // if gt.in_this_processor
		}// I1
	}// T1


	//------------
	// for test
	//------------
/*	cout << setprecision(3);
	for(int i=0; i<nnrg_now; i++)

	for(int ik=0; ik<kv.nkstot; ++ik)
	{
		for(int ib=0; ib<NBANDS; ++ib)
		{
			cout << " ik=" << ik << " ib=" << ib << " occ=" << wf.wg(ik,ib) << " e=" << wf.ekb[ik][ib] << endl;
		}
	}

	for(int i=0; i<10; i++)
	{
		if(DM_R[0][i]>1.0e-8)
		{
			cout << " i=" << i << " DM_R=" << DM_R[0][i] << endl;
		}
	}
*/
	for(int i=0; i<NSPIN; ++i)
		delete[] DM_ATOM[i];
	delete[] DM_ATOM;
	delete[] WFC_PHASE;
	RA.delete_grid();//xiaohui add 2015-02-04
	timer::tick("LCAO_Charge","cal_dk_k",'F');	
	return;
}

//-------------------------------------------------------------
// NOTE:
// How to improve the speed of calculation of density matrix.
// There are two ways to calculate density matrix.
// 1) After diagonalization, we get c_nk, where n is band index
// and k in k index.
// we distribute c_n_mu (mu is orbital) to different processors 
// to generate density matrix by the formula:
// rhodm_mu_nu=\sum_{n}c_n_mu * c_n_nu
// 2) After diagonalization, we generate density matrix first
// rhodm_mu_nu=\sum_{n}c_n_mu * c_n_nu and then we distribute
// density matrix.
//
// I am not sure which one is slower, but I guess by now
// the first one we use seems to be slow for large system
// because we need the data to distribute is increased
// with O(N^2) (first N from n, bands, second N from mu,
// orbitals), so that's not a good method.
//
// Another advantage we haven't taken is the density matrix
// is symmetried.
//
// 2014-05-18 Mohan
//
//--------------------------------------------------------------
void Local_Orbital_Charge::cal_dk_gamma(void)
{
	 TITLE("Local_Orbital_Charge","cal_density_kernal");
	timer::tick("LocalOrbital_Charge","cal_dk_gamma",'F');

	assert(NSPIN==kv.nks);

	//xiaohui add 2014-06-18
	//for(int is=0; is<NSPIN; is++)
	//{
	//	for (int ib=0; ib<NBANDS; ib++)
	//	{
	//		const double w1 = wf.wg(is, ib);
	//		ofs_running <<"w1[" <<is <<"]["<<ib <<"]: "<< w1 << endl; 
	//	}
	//}

	if(BFIELD)
	{
		for(int is=0; is<NSPIN; is++)
		{
			for (int i=0; i<lgd_now; i++)
			{
				ZEROS(this->DM_B[is][i], lgd_now);
			}
		}
		for(int is=0; is<NSPIN; is++)
		{
			for (int i=0; i<NLOCAL; i++)
			{
				const int mu_local = GridT.trace_lo[i];
				if ( mu_local >= 0)
				{
					// set a pointer.
					complex<double> *alpha = this->DM_B[is][mu_local];
					for (int j=i; j<NLOCAL; j++)
					{
						const int nu_local = GridT.trace_lo[j];
						if ( nu_local >= 0)
						{
							for (int ib=0; ib<NBANDS; ib++)
							{
								const double w1 = wf.wg(is, ib);
								if(w1>0)
								{
									// dm = \sum ( wg * c[ib][mu] * c[ib][nu] )
									// dm saved accordint to sub-FFT box.
									alpha[nu_local] += w1 * (conj(LOWF.WFC_GAMMA_B[is][ib][mu_local]) * LOWF.WFC_GAMMA_B[is][ib][nu_local]).real();
								}
							}
						}
					}
				}
			}
		}
	}
#ifdef __MPI //2015-09-06, xiaohui
	else
	{
		for(int is=0; is<NSPIN; is++)
		{
			for (int i=0; i<lgd_now; i++)
			{
				ZEROS(this->DM[is][i], lgd_now);
			}
		}
		//for(int is=0; is<NSPIN; is++)
		//{
		//	for (int i=0; i<NLOCAL; i++)
		//	{
		//		ZEROS(this->DM[is][i], NLOCAL);
		//	}
		//}

		//double **out_DM;
		//out_DM = new double *[NSPIN];
		//for(int is=0; is<NSPIN; is++)
		//{
		//	out_DM[is] = new double[NLOCAL*NLOCAL];
		//	ZEROS(out_DM[is], NLOCAL*NLOCAL);
		//}

		int nprocs,myid;
		MPI_Status status;
		MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
		MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

		int local_band=NBANDS/DSIZE;
	 		if (DRANK<NBANDS%DSIZE) local_band=local_band+1;

		int *local_bands;
	 		local_bands = new int[DSIZE];
		ZEROS(local_bands, DSIZE);

	 		int lastband_in_proc = 0;
	 		int lastband_number = 0;
	 		int count_bands = 0;
		for (int i=0; i<DSIZE; i++)
	 		{
		  		if (i<NBANDS%DSIZE)
		  		{
							local_bands[i]=NBANDS/DSIZE+1;
		  		}
		  		else
		  		{
							local_bands[i]=NBANDS/DSIZE;
		  		}
		  		count_bands += local_bands[i];
		  		if (count_bands >= NBANDS)
		  		{
							lastband_in_proc = i;
							lastband_number = NBANDS - (count_bands - local_bands[i]);
							break;
		  		}
	 		}
		//double *w1 ;
		//w1 = new double[local_band] ;
		//ZEROS(w1, local_band);
		double** w1;
		w1 = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			w1[is] = new double[local_band];
			ZEROS(w1[is], local_band);
		}

		int total_bands = 0;

		time_t time_wg2w1_start = time(NULL);
		double Time_Wg2w1_Start = (double)clock();

		int Total_Bands = 0;
		int Nbands = -1;
		for(int i=0; i <= lastband_in_proc; i++)
		{
			Nbands = local_bands[i];
			if(myid == i)
			{
				for(int is=0; is<NSPIN; is++)
				{
					for (int ib=0; ib< Nbands; ib++)
					{
						w1[is][ib] = wf.wg(is, ib + Total_Bands);
					}
				}
			}
			Total_Bands += Nbands;
		}

		//double* Z_wg;
		//Z_wg = new double[local_band * NLOCAL];
		//ZEROS(Z_wg, local_band * NLOCAL);
		double** Z_wg;
		Z_wg = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			Z_wg[is] = new double[local_band * NLOCAL];
			ZEROS(Z_wg[is], local_band * NLOCAL);
		}

		time_t time_Z_wg_start = time(NULL);
		double Time_Z_Wg_Start = (double)clock();

		//stringstream ss;
		//ss<<"wavefunction_on_processor_"<<MY_RANK+1;
		////ss<<"FULL_DM";
		//ofstream ofs;
		//ofs.open(ss.str().c_str());

		//stringstream ss1;
		//ss1<<"weight_on_processor_"<<MY_RANK+1;
		////ss<<"FULL_DM";
		//ofstream ofs1;
		//ofs1.open(ss1.str().c_str());

		if(myid <= lastband_in_proc)
		{
			for(int is=0; is<NSPIN; is++)
			{
				for(int i=0; i<local_bands[myid]; i++)
				{
					//if(i != 0 && i%10 == 0) ofs1<<endl;
					//ofs1<<w1[is][i]<<"\t";
					for(int n=0; n<NLOCAL; n++)
					{
						Z_wg[is][i + n*local_bands[myid]] = ParaO.Z_LOC[is][i + n*local_bands[myid]]* w1[is][i];
						//if(n != 0 && n%10 == 0) ofs<<endl;
						//ofs<<ParaO.Z_LOC[is][i + n*local_bands[myid]]<<"\t";
					}
				}
			}
		}
//cout<<"Z_wg[0]: "<<ParaO.Z_LOC[0]<<"\t"<<ParaO.Z_LOC[1]<<"\t"<<ParaO.Z_LOC[2]<<endl;
//cout<<"Z_wg[0]: "<<Z_wg[0]<<"\t"<<Z_wg[1]<<"\t"<<Z_wg[2]<<endl;

		//delete[] w1;

		//stringstream ss;
		////ss<<"EDM_"<<MY_RANK+1;
		//ss<<"FULL_DM";
		//ofstream ofs;
		//ofs.open(ss.str().c_str());
		////ofs<<"local_band: "<<local_band<<"\t"<<"local_bands[id]: "<<local_bands[myid]<<endl;

		for(int is=0; is<NSPIN; is++)
		{
			delete[] w1[is];
		}
		delete[] w1;
		delete[] local_bands;

		time_t time_Z_wg_end = time(NULL);
		double Time_Z_Wg_End =(double)clock();
		double time_Z_wg = difftime(time_Z_wg_start,time_Z_wg_end);
		//OUT_TIME("calculate Z_wg",time_wg_start,time_wg_end);

		//int local_band = local_band;
		double coll_time = 0.0;
		double Time_Cal_Start = (double)clock();

		int row_col = NLOCAL/300;
		if(row_col >= 1)
		{
			//double *ZW_300;
			//ZW_300 = new double[300*local_band];
			////ZEROS(ZW_300, 300*local_band);
			double** ZW_300;
			ZW_300 = new double*[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				ZW_300[is] = new double[300 * local_band];
				ZEROS(ZW_300[is], 300 * local_band);
			}
			//double *Z_300;
			//Z_300 = new double[300*local_band];
			////ZEROS(Z_300, 300*local_band);
			double** Z_300;
			Z_300 = new double*[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				Z_300[is] = new double[300 * local_band];
				ZEROS(Z_300[is], 300 * local_band);
			}

			if(NLOCAL%300 != 0)	row_col += 1;
			int row_global = 0;
			int row_index = 0;
			int col_global = 0;
			int col_index = 0;
			for(int row_count=1; row_count<=row_col; row_count++)
			{
				row_global = row_count*300;
				if(row_global <= NLOCAL)
				{
					for(int is=0; is<NSPIN; is++)
					{
						//ZEROS(ZW_300, 300*local_band);
						for(int i_row=0; i_row<300; i_row++)
						{
							row_index = i_row +(row_count-1)*300;
							for(int ib=0; ib<local_band; ib++)
							{
								ZW_300[is][ib + i_row*local_band] = Z_wg[is][ib + row_index*local_band] ;
							}
						}
					}
					for(int col_count=1; col_count<=row_col; col_count++)
					{
						col_global = col_count*300;
						if(col_global <= NLOCAL)
						{
							//ZEROS(Z_300, 300*local_band);
							double **rho_300;
							rho_300 = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_300[is] = new double[300*300];
								ZEROS(rho_300[is], 300*300);

								for(int i_col=0; i_col<300; i_col++)
								{
									col_index = i_col +(col_count-1)*300;
									for(int ib=0; ib<local_band; ib++)
									{
										Z_300[is][ib + i_col*local_band] = ParaO.Z_LOC[is][ib + col_index*local_band] ;
									}
								}

								//for(int i_row=0; i_row<300; i_row++)
								//{
								//	for(int i_col=0; i_col<300; i_col++)
								//	{
								//		for(int ib=0; ib<local_band; ib++)
								//		{
								//			rho_300[is][i_col + i_row*300] += ZW_300[is][ib + i_row*local_band]  * Z_300[is][ib + i_col*local_band];
								//		}
								//	}
								//}
					 						double alpha=1;
					 						double beta=0;
					 						char transa = 'T';
					 						char transb = 'N';
					 						int M = 300;
					 						int N = 300;
					 						//int K = local_band;
					 						////int lda = local_band;
					 						////int ldb = local_band;
					 						int K = local_band;
					 						int lda = local_band;
					 						int ldb = local_band;
					 						int ldc = 300;

								dgemm_(&transa, &transb, &M, &N, &K,
				  								&alpha, &ZW_300[is][0], &lda, &Z_300[is][0], &ldb,
			 									&beta, &rho_300[is][0], &ldc);
								double Time_Coll_Start = (double)clock();
								//MPI_Barrier(MPI_COMM_WORLD);
								MPI_Barrier(DIAG_HPSEPS_WORLD);
								Parallel_Reduce::reduce_double_all( rho_300[is], 300*300);
								double Time_Coll_End = (double)clock();

								//coll_time += (Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC);
								coll_time += (Time_Coll_End-Time_Coll_Start);
								if(GAMMA_ONLY_LOCAL)
		  							{
									//for(int is=0; is<NSPIN; is++)
									//{
									for(int i_row=0; i_row<300; i_row++)
									{
										row_index = i_row +(row_count-1)*300;
										int row_mu = GridT.trace_lo[row_index];
										//if(row_mu >= 0)
										//{
											for(int i_col=0; i_col<300; i_col++)
											{
												col_index = i_col +(col_count-1)*300;
												int col_nu = GridT.trace_lo[col_index];
//out_DM[is][col_index + row_index*NLOCAL] = rho_300[is][i_row + i_col*300];
												if(row_mu >= 0 && col_nu >= 0)
												//if(col_nu >= 0)
												{
													//this->DM[is][row_mu][col_nu] = rho_300[i_col + i_row*300];
													//this->DM[is][row_mu][col_nu] = rho_300[is][i_col + i_row*300];
													this->DM[is][row_mu][col_nu] = rho_300[is][i_row + i_col*300];
												}
											}
										//}
									}
								}
								delete[] rho_300[is];
							}
							delete[] rho_300;
						}
						else
						{
							int col_remain = 0;
							int col_index2 = 0;
							col_remain = NLOCAL - (col_count - 1)*300;
							//double *Z_remain = new double[col_remain * local_band];
							//ZEROS(Z_remain, col_remain*local_band);
							double** Z_remain;
							Z_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								Z_remain[is] = new double[col_remain * local_band];
								ZEROS(Z_remain[is], col_remain * local_band);
							}

							//double *rho_300_remain = new double[300*col_remain];
							//ZEROS(rho_300_remain, 300*col_remain);
							double** rho_300_remain;
							rho_300_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_300_remain[is] = new double[300*col_remain];
								ZEROS(rho_300_remain[is], 300*col_remain);

								for(int i_col=0; i_col<col_remain; i_col++)
								{
									col_index2 = i_col +(col_count-1)*300;
									for(int ib=0; ib<local_band; ib++)
									{
										Z_remain[is][ib + i_col*local_band] = ParaO.Z_LOC[is][ib + col_index2*local_band] ;
									}
								}

								//xiaohui modify 2015-03-25
								//for(int i_row=0; i_row<300; i_row++)
								//{
								//	for(int i_col=0; i_col<col_remain; i_col++)
								//	{
								//		//for(int ib=0; ib<local_band; ib++)
								//		for(int ib=0; ib<local_band; ib++)
								//		{
								//			//rho_300_remain[i_col + i_row*col_remain] += ZW_300[ib + i_row*local_band]  * Z_remain[ib + i_col*local_band];
								//			rho_300_remain[is][i_col + i_row*col_remain] += ZW_300[is][ib + i_row*local_band]  * Z_remain[is][ib + i_col*local_band];
								//		}
								//	}
								//}
					 			double alpha=1;
					 			double beta=0;
					 			char transa = 'T';
					 			char transb = 'N';
					 			int M = 300;
					 			int N = col_remain;
					 			////int K = local_band;
					 			////int lda = local_band;
					 			////int ldb = local_band;
					 			int K = local_band;
					 			int lda = local_band;
					 			int ldb = local_band;
					 			int ldc = 300;
					 			//int ldc = col_remain;

								dgemm_(&transa, &transb, &M, &N, &K,
				  					&alpha, &ZW_300[is][0], &lda, &Z_remain[is][0], &ldb,
			 						&beta, &rho_300_remain[is][0], &ldc);
								double Time_Coll_Start = (double)clock();
								//MPI_Barrier(MPI_COMM_WORLD);
								MPI_Barrier(DIAG_HPSEPS_WORLD);
								Parallel_Reduce::reduce_double_all( rho_300_remain[is], 300*col_remain);
								double Time_Coll_End = (double)clock();

								//coll_time += (Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC);
								coll_time += (Time_Coll_End-Time_Coll_Start);
//cout<<"300*col_remain time: "<<(Time_Coll_End-Time_Coll_Start)/(C
								//MPI_Barrier(MPI_COMM_WORLD);
								//Parallel_Reduce::reduce_double_all( rho_300_remain, 300*col_remain);
								if(GAMMA_ONLY_LOCAL)
		  							{
									//for(int is=0; is<NSPIN; is++)
									//{
									for(int i_row=0; i_row<300; i_row++)
									{
										row_index = i_row +(row_count-1)*300;
										int row_mu = GridT.trace_lo[row_index];
										//if(row_mu >= 0)
										//{
											for(int i_col=0; i_col<col_remain; i_col++)
											{
												col_index = i_col +(col_count-1)*300;
												int col_nu = GridT.trace_lo[col_index];
//out_DM[is][col_index + row_index*NLOCAL] = rho_300_remain[is][i_row + i_col*300];
												if(row_mu >= 0 && col_nu >= 0)
												//if(col_nu >= 0)
												{
													//this->DM[is][row_mu][col_nu] = rho_300_remain[i_col + i_row*col_remain];
													//this->DM[is][row_mu][col_nu] = rho_300_remain[is][i_col + i_row*col_remain];
													this->DM[is][row_mu][col_nu] = rho_300_remain[is][i_row + i_col*300];
												}
											}
										//}
									}
								}
								delete[] Z_remain[is];
								delete[] rho_300_remain[is];
							}
							delete[] Z_remain;
							delete[] rho_300_remain;
						}
					}
				}
				else
				{
					int row_remain = 0;
					int row_index2 = 0;
					row_remain = NLOCAL - (row_count - 1)*300;
					//double *Z_row_remain = new double[row_remain * local_band];
					//ZEROS(Z_row_remain, row_remain*local_band);
					double** Z_row_remain;
					Z_row_remain = new double*[NSPIN];
					for(int is=0; is<NSPIN; is++)
					{
						Z_row_remain[is] = new double[row_remain * local_band];
						ZEROS(Z_row_remain[is], row_remain * local_band);
					}

					//double *rho_remain_300 = new double[row_remain*300];
					//ZEROS(rho_remain_300, row_remain*300);
					//double** rho_remain_300;
					//rho_remain_300 = new double*[NSPIN];
					for(int is=0; is<NSPIN; is++)
					{
						//rho_remain_300[is] = new double[row_remain * 300];
						//ZEROS(rho_remain_300[is], row_remain * 300);

						for(int i_row=0; i_row<row_remain; i_row++)
						{
							row_index2 = i_row +(row_count-1)*300;
							for(int ib=0; ib<local_band; ib++)
							{
								Z_row_remain[is][ib + i_row*local_band] = Z_wg[is][ib + row_index2*local_band] ;
							}
						}
					}

					for(int col_count=1; col_count<=row_col; col_count++)
					{
						col_global = col_count*300;
						if(col_global <= NLOCAL)
						{
							double** rho_remain_300;
							rho_remain_300 = new double*[NSPIN];

							for(int is=0; is<NSPIN; is++)
							{
								rho_remain_300[is] = new double[row_remain * 300];
								ZEROS(rho_remain_300[is], row_remain * 300);

								//ZEROS(Z_300, 300*local_band);
								for(int i_col=0; i_col<300; i_col++)
								{
									col_index = i_col +(col_count-1)*300;
									for(int ib=0; ib<local_band; ib++)
									{
										Z_300[is][ib + i_col*local_band] = ParaO.Z_LOC[is][ib + col_index*local_band] ;
									}
								}

								//xiaohui modify 2015-03-25
								//for(int i_row=0; i_row<row_remain; i_row++)
								//{
								//	for(int i_col=0; i_col<300; i_col++)
								//	{
								//		//for(int ib=0; ib<local_band; ib++)
								//		for(int ib=0; ib<local_band; ib++)
								//		{
								//			rho_remain_300[is][i_col + i_row*300] += Z_row_remain[is][ib + i_row*local_band]  * Z_300[is][ib + i_col*local_band];
								//		}
								//	}
								//}
					 			double alpha=1;
					 			double beta=0;
					 			char transa = 'T';
					 			char transb = 'N';
					 			int M = row_remain;
					 			////int M = 300;
					 			int N = 300;
					 			////int K = local_band;
					 			////int lda = local_band;
					 			////int ldb = local_band;
					 			int K = local_band;
					 			int lda = local_band;
					 			int ldb = local_band;
					 			////int ldc = 300;
					 			int ldc = row_remain;

								dgemm_(&transa, &transb, &M, &N, &K,
				  					&alpha, &Z_row_remain[is][0], &lda, &Z_300[is][0], &ldb,
			 						&beta, &rho_remain_300[is][0], &ldc);
								double Time_Coll_Start = (double)clock();
								//MPI_Barrier(MPI_COMM_WORLD);
								MPI_Barrier(DIAG_HPSEPS_WORLD);
								Parallel_Reduce::reduce_double_all( rho_remain_300[is], row_remain*300);
								double Time_Coll_End = (double)clock();

//cout<<"row_remain*300 time: "<<(Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC)<<endl;
								//coll_time += (Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC);
								coll_time += (Time_Coll_End-Time_Coll_Start);
								//MPI_Barrier(MPI_COMM_WORLD);
								//Parallel_Reduce::reduce_double_all( rho_300_remain2, row_remain*300);
								if(GAMMA_ONLY_LOCAL)
		  							{
									//for(int is=0; is<NSPIN; is++)
									//{
										for(int i_row=0; i_row<row_remain; i_row++)
										{
											row_index = i_row +(row_count-1)*300;
											int row_mu = GridT.trace_lo[row_index];
											//if(row_mu >= 0)
											//{
												for(int i_col=0; i_col<300; i_col++)
												{
													col_index = i_col +(col_count-1)*300;
													int col_nu = GridT.trace_lo[col_index];
//out_DM[is][col_index + row_index*NLOCAL] = rho_remain_300[is][i_row + i_col*row_remain];
													if(row_mu >= 0 && col_nu >= 0)
													//if(col_nu >= 0)
													{
														//this->DM[is][row_mu][col_nu] = rho_300_remain2[i_col + i_row*300];
														//this->DM[is][row_mu][col_nu] = rho_remain_300[is][i_col + i_row*300];
														this->DM[is][row_mu][col_nu] = rho_remain_300[is][i_row + i_col*row_remain];
													}
												}
											//}
										}
								}
							}
							for(int is=0; is<NSPIN; is++)
							{
								delete[] rho_remain_300[is];
							}
							delete[] rho_remain_300;
						}
						else
						{
							int col_remain = 0;
							int col_index2 = 0;
							col_remain = NLOCAL - (col_count - 1)*300;
							//double *Z_col_remain = new double[col_remain * local_band];
							//ZEROS(Z_col_remain, col_remain*local_band);
							double** Z_col_remain;
							Z_col_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								Z_col_remain[is] = new double[col_remain * local_band];
								ZEROS(Z_col_remain[is], col_remain * local_band);
							}
							//double *rho_remain_remain = new double[row_remain*col_remain];
							//ZEROS(rho_remain_remain, row_remain*col_remain);
							double** rho_remain_remain;
							rho_remain_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_remain_remain[is] = new double[row_remain * col_remain];
								ZEROS(rho_remain_remain[is], row_remain * col_remain);
							//2014-09-22 12:04
								for(int i_col=0; i_col<col_remain; i_col++)
								{
									col_index2 = i_col +(col_count-1)*300;
									for(int ib=0; ib<local_band; ib++)
									{
										Z_col_remain[is][ib + i_col*local_band] = ParaO.Z_LOC[is][ib + col_index2*local_band] ;
									}
								}

								//for(int i_row=0; i_row<row_remain; i_row++)
								//{
								//	for(int i_col=0; i_col<col_remain; i_col++)
								//	{
										//for(int ib=0; ib<local_band; ib++)
								//		for(int ib=0; ib<local_band; ib++)
								//		{
								//			rho_remain_remain[is][i_col + i_row*col_remain] += Z_row_remain[is][ib + i_row*local_band]  * Z_col_remain[is][ib + i_col*local_band];
								//		}
								//	}
								//}
					 			double alpha=1;
					 			double beta=0;
					 			char transa = 'T';
					 			char transb = 'N';
					 			////int M = 300;
					 			int M = row_remain;
					 			int N = col_remain;
					 			////int K = local_band;
					 			////int lda = local_band;
					 			////int ldb = local_band;
					 			int K = local_band;
					 			int lda = local_band;
					 			int ldb = local_band;
					 			////int ldc = 300;
					 			int ldc = row_remain;

								dgemm_(&transa, &transb, &M, &N, &K,
				  					&alpha, &Z_row_remain[is][0], &lda, &Z_col_remain[is][0], &ldb,
			 						&beta, &rho_remain_remain[is][0], &ldc);
								double Time_Coll_Start = (double)clock();
								//MPI_Barrier(MPI_COMM_WORLD);
								MPI_Barrier(DIAG_HPSEPS_WORLD);
								Parallel_Reduce::reduce_double_all( rho_remain_remain[is], row_remain*col_remain);
								double Time_Coll_End = (double)clock();
//cout<<"row_remain*col_remain time: "<<(Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC)<<endl;

								//coll_time += (Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC);
								coll_time += (Time_Coll_End-Time_Coll_Start);
								//MPI_Barrier(MPI_COMM_WORLD);
								//Parallel_Reduce::reduce_double_all( rho_remain_remain, row_remain*col_remain);
								if(GAMMA_ONLY_LOCAL)
		  							{
									//for(int is=0; is<NSPIN; is++)
									//{
										for(int i_row=0; i_row<row_remain; i_row++)
										{
											row_index = i_row +(row_count-1)*300;
											int row_mu = GridT.trace_lo[row_index];
											//if(row_mu >= 0)
											//{
												for(int i_col=0; i_col<col_remain; i_col++)
												{
													col_index = i_col +(col_count-1)*300;
													int col_nu = GridT.trace_lo[col_index];
//out_DM[is][col_index + row_index*NLOCAL] = rho_remain_remain[is][i_row + i_col*row_remain];
													if(row_mu >= 0 && col_nu >= 0)
													//if(col_nu >= 0)
													{
														//this->DM[is][row_mu][col_nu] = rho_remain_remain[i_col + i_row*col_remain];
														//this->DM[is][row_mu][col_nu] = rho_remain_remain[is][i_col + i_row*col_remain];
														this->DM[is][row_mu][col_nu] = rho_remain_remain[is][i_row + i_col*row_remain];
													}
												}
											//}
										}
								}
								delete[] Z_col_remain[is];
								delete[] rho_remain_remain[is];
							}
							delete[] Z_col_remain;
							delete[] rho_remain_remain;
							
						}
					}
					for(int is=0; is<NSPIN; is++)
					{
						delete[] Z_row_remain[is];
					}
					delete[] Z_row_remain;
				}
			}
			for(int is=0; is<NSPIN; is++)
			{
				delete[] ZW_300[is];
				delete[] Z_300[is];
			}
			delete[] ZW_300;
			delete[] Z_300;
		}
		else
		{
			double** rho_NLOCAL_NLOCAL;
			rho_NLOCAL_NLOCAL = new double*[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				rho_NLOCAL_NLOCAL[is] = new double[NLOCAL*NLOCAL];
				ZEROS(rho_NLOCAL_NLOCAL[is], NLOCAL*NLOCAL);
			}
			for(int is=0; is<NSPIN; is++)
			{
				for(int i_row=0; i_row<NLOCAL; i_row++)
				{
					for(int i_col=0; i_col<NLOCAL; i_col++)
					{
						//for(int ib=0; ib<local_band; ib++)
						for(int ib=0; ib<local_band; ib++)
						{
							//rho_300_remain[i_col + i_row*col_remain] += ZW_300[ib + i_row*local_band]  * Z_remain[ib + i_col*local_band];
							rho_NLOCAL_NLOCAL[is][i_col + i_row*NLOCAL] += Z_wg[is][ib + i_row*local_band]  * ParaO.Z_LOC[is][ib + i_col*local_band];
						}
					}
				}
				double Time_Coll_Start = (double)clock();
				MPI_Barrier(DIAG_HPSEPS_WORLD);
				Parallel_Reduce::reduce_double_all( rho_NLOCAL_NLOCAL[is], NLOCAL*NLOCAL);
				double Time_Coll_End = (double)clock();
				//coll_time += (Time_Coll_End-Time_Coll_Start)/(CLOCKS_PER_SEC);
				coll_time += (Time_Coll_End-Time_Coll_Start);
				//MPI_Barrier(MPI_COMM_WORLD);
				//Parallel_Reduce::reduce_double_all( rho_remain_remain, row_remain*col_remain);
				if(GAMMA_ONLY_LOCAL)
		  			{
					//for(int is=0; is<NSPIN; is++)
					//{
					for(int i_row=0; i_row<NLOCAL; i_row++)
					{
						int row_mu = GridT.trace_lo[i_row];
						//if(row_mu >= 0)
						//{
							for(int i_col=0; i_col<NLOCAL; i_col++)
							{
								int col_nu = GridT.trace_lo[i_col];
//out_DM[is][i_col + i_row*NLOCAL] = rho_NLOCAL_NLOCAL[is][i_col + i_row*NLOCAL];
								if(row_mu >= 0 && col_nu >= 0)
								//if(col_nu >= 0)
								{
									this->DM[is][row_mu][col_nu] = rho_NLOCAL_NLOCAL[is][i_col + i_row*NLOCAL];
								}
							}
						//}
					}
				}
				delete[] rho_NLOCAL_NLOCAL[is];
			}
			delete[] rho_NLOCAL_NLOCAL;
		}
		//stringstream ss;
		////ss<<"EDM_"<<MY_RANK+1;
		//ss<<"FULL_DM";
		//ofstream ofs;
		//ofs.open(ss.str().c_str());
		////ofs<<"ParaO.nrow: "<<ParaO.nrow<<"\t"<<"ParaO.ncol: "<<ParaO.ncol<<endl;

		//for(int is=0; is<NSPIN; is++)
		//{
		//	for(int i=0; i<NLOCAL*NLOCAL; i++)
		//	{
		//		if(i != 0 && i%10 == 0) ofs<<endl;
		//		ofs<<out_DM[is][i]<<"\t";
		//	}
		//}

		for(int is=0; is<NSPIN; is++)
		{
			delete[] Z_wg[is];
			//delete[] out_DM[is];
		}
		delete[] Z_wg;
		//delete[] out_DM;
	}
#endif //2015-09-06, xiaohui
#ifndef __MPI //2015-09-06, xiaohui
	else
	{
		for(int is=0; is<NSPIN; is++)
		{
			for (int i=0; i<lgd_now; i++)
			{
				ZEROS(this->DM[is][i], lgd_now);
			}
		}
		for(int is=0; is<NSPIN; is++)
		{
			for (int i=0; i<NLOCAL; i++)
			{
				const int mu_local = GridT.trace_lo[i];
//				ofs_running << " mu_local=" << mu_local << endl;
				if ( mu_local >= 0)
				{
					// set a pointer.
					double *alpha = this->DM[is][mu_local];
					for (int j=i; j<NLOCAL; j++)
					{
						const int nu_local = GridT.trace_lo[j];
						if ( nu_local >= 0)
						{
							for (int ib=0; ib<NBANDS; ib++)
							{
								const double w1 = wf.wg(is, ib);
								//ofs_running << " w1=" << w1 << endl;
								if(w1>0)
								{
									// dm = \sum ( wg * c[ib][mu] * c[ib][nu] )
									// dm saved accordint to sub-FFT box.
									alpha[nu_local] += w1 * LOWF.WFC_GAMMA[is][ib][mu_local] 
									* LOWF.WFC_GAMMA[is][ib][nu_local];

								}//w1
							}//ib
						}//nu_local
					}//j
				}//mu_local
			}//i
		}//is
	}
#endif //2015-09-06, xiaohui

	//----------
	// for test
	//----------
	/*
	ofs_running << " WFC2 " << endl;
	for(int i=0; i<NBANDS; ++i)
	{
		ofs_running << " wf.eg=" << wf.wg(0,i) << endl;
		for(int j=0; j<NLOCAL; ++j)
		{
			const int mu = GridT.trace_lo[j];
			if(mu>=0)
			{
				if( abs(LOWF.WFC_GAMMA[0][i][mu] > 1.0e-8) )
				ofs_running << setw(5) << i+1 << setw(8) << j+1 << setw(15) << LOWF.WFC_GAMMA[0][i][mu] << endl; 
			}
		}
	}
	*/


	// @@@@@@@
	// test
	// @@@@@@@
	/*
	ofs_running << " LOWF.DM" << endl;
	for(int i=0; i<NBANDS; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			if( GridT.trace_lo[j] >= 0)
			{
//				if( abs(DM[0][i][j]) > 1.0e-8 )
				{
					ofs_running << setw(5) << i+1 << setw(8) << j+1 << setw(15) << this->DM[0][i][GridT.trace_lo[j]] << endl;
				}
			}
		}
	}
	*/

	timer::tick("LocalOrbital_Charge","cal_dk_gamma",'F');
	 return;
}



//-------------------------------------------------
// NOTE for Local_Orbital_Charge::write_dm
// I will give an example here, suppose we have a 4*4 
// density matrix (symmetry) which is
// 1.1 2.3 3.6 4.2
// 2.3 5.2 7.1 8.9
// 3.6 7.1 1.3 3.2  
// 4.2 8.9 3.2 2.4
// we use two processors, each one has 3 orbitals
// processor 1 has orbital index 1, 2, 4:
// ('no' means no value on this processor)
// 1.1 2.3 no  4.2
// 2.3 5.2 no  8.9
// no  no  no  no	
// 4.2 8.9 no  2.4
// processor 2 has orbital index 1, 3, 4;
// 1.1 no  3.6 4.2
// no  no  no  no 
// 3.6 no  1.3 3.2  
// 4.2 no  3.2 2.4
// now we want to reduce them and print out,
// we plan to reduce them one row by one row,
// then for the first row, we need to set the
// temparary array to 4 (NLOCAL in code),
// then we reduce first row, it will become
// 2.2 2.3 3.6 8.4,
// we notice that the first element and fourth
// element is doubled, that's because the density
// may have overlap, so we need to first count
// for each element, how many times it is duplicated
// on other processors, this is why there is
// a 'count' integer array in the code.
// UPDATED BY MOHAN 2014-05-18
void Local_Orbital_Charge::write_dm(const int &is, const int &iter, const string &fn, const int &precision)
{
	TITLE("Local_Orbital_Charge","write_dm");

	 if (out_dm==0)
	 {
		  return;
	 }
	 else if(iter % out_dm != 0)
	 {
		  return; 
	 }
	timer::tick("Local_Orbital_Charge","write_dm");

	 time_t start, end;
	 ofstream ofs;

	 if(MY_RANK==0)
	 {
		  start = time(NULL);

		  ofs.open(fn.c_str());
		  if (!ofs)
		  {
				WARNING("Charge::write_rho","Can't create Charge File!");
		  }

		  //ofs_running << "\n Output charge file." << endl;

		  ofs << ucell.latName << endl;//1
		  ofs << " " << ucell.lat0 * BOHR_TO_A << endl;
		  ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		  ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		  ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;
		  for(int it=0; it<ucell.ntype; it++)
		  {
				ofs << " " << ucell.atoms[it].label;
		  }
		  ofs << endl;
		  for(int it=0; it<ucell.ntype; it++)
		  {
				ofs << " " << ucell.atoms[it].na;
		  }
		  ofs << endl;
		  ofs << "Direct" << endl;

		  for(int it=0; it<ucell.ntype; it++)
		  {
			Atom* atom = &ucell.atoms[it];
			ofs << setprecision(15);
				for(int ia=0; ia<ucell.atoms[it].na; ia++)
				{
					 ofs << " " << atom->taud[ia].x
						  << " " << atom->taud[ia].y
						  << " " << atom->taud[ia].z << endl;
				}
		  }

		ofs << "\n " << NSPIN;
		if(NSPIN==1)
		{
			ofs << "\n " << en.ef << " (fermi energy)";
		}
		else if(NSPIN==2)
		{
			if(is==0)ofs << "\n " << en.ef_up << " (fermi energy for spin=1)";
			else if(is==1)ofs << "\n " << en.ef_dw << " (fermi energy for spin=2)";
		}
		else
		{
			WARNING_QUIT("write_rho","check nspin!");
		}

	
		ofs << "\n  " << NLOCAL << " " << NLOCAL << endl;

		  ofs << setprecision(precision);
		  ofs << scientific;

	 }

	//ofs << "\n " << GAMMA_ONLY_LOCAL << " (GAMMA ONLY LOCAL)" << endl;
#ifndef __MPI
	if(GAMMA_ONLY_LOCAL)
	{
		for(int i=0; i<NLOCAL; ++i)
		{
			for(int j=0; j<NLOCAL; ++j)
			{
				if(j%8==0) ofs << "\n";
				ofs << " " << this->DM[is][i][j];
			}
		}
	}
	else
	{
		WARNING_QUIT("write_dm","not ready yet");
		ofs << " " << LNNR.nnrg << " (nnrg)" << endl;
		for(int i=0; i<LNNR.nnrg; ++i)
		{
			if(i%8==0) ofs << "\n";
			ofs << " " << this->DM_R[is][i];
		}
	}
#else
	if(GAMMA_ONLY_LOCAL)
	{
		//xiaohui modify 2014-06-18
		
		double* tmp = new double[NLOCAL];
		int* count = new int[NLOCAL];
		for (int i=0; i<NLOCAL; ++i)
		{
			// when reduce, there may be 'redundance', we need to count them.
			ZEROS(count, NLOCAL);
			const int mu = GridT.trace_lo[i];
			if (mu >= 0)
			{
				for (int j=0; j<NLOCAL; ++j)
				{
					const int nu = GridT.trace_lo[j];
					if (nu >= 0)
					{
						count[j]=1;	
					}
				}
			}
			Parallel_Reduce::reduce_int_all( count, NLOCAL );

			// reduce the density matrix for 'i' line.
			ZEROS(tmp, NLOCAL);
			if (mu >= 0)
			{
				for (int j=0; j<NLOCAL; j++)
				{
					const int nu = GridT.trace_lo[j];
					if (nu >=0)
					{
						tmp[j] = DM[is][mu][nu];
						//ofs_running << " dmi=" << i << " j=" << j << " " << DM[is][mu][nu] << endl;
					}
				}
			}
			Parallel_Reduce::reduce_double_all( tmp, NLOCAL );

			if(MY_RANK==0)
			{
				for (int j=0; j<NLOCAL; j++)
				{
					if(j%8==0) ofs << "\n";
					if(count[j]>0)
					{
						ofs << " " << tmp[j]/(double)count[j];
					}
					else
					{
						ofs << " 0"; 
					}
				}
			}
		}
		delete[] tmp;
		delete[] count;
		
		//xiaohui add 2014-06-18
		//for(int i=0; i<NLOCAL; ++i)
		//{
		//	for(int j=0; j<NLOCAL; ++j)
		//	{
		//		if(j%8==0) ofs << "\n";
		//		ofs << " " << this->DM[is][i][j];
		//	}
		//}

	}
	else
	{
		ofs << " " << LNNR.nnrg << " (nnrg)" << endl;
		WARNING_QUIT("local_orbital_charge","not ready to output DM_R");
	}
#endif
	 if(MY_RANK==0)
	 {
		  end = time(NULL);
		  OUT_TIME("write_rho",start,end);
		  ofs.close();
	 }
	timer::tick("Local_Orbital_Charge","write_dm");

	return;
}


void Local_Orbital_Charge::read_dm(const int &is, const string &fn)
{
	TITLE("Local_Orbital_Charge","read_dm");
	timer::tick("Local_Orbital_Charge","read_dm");

	ofs_running << "\n processor 0 is reading density matrix from file < " << fn << " > " << endl;
	//xiaohui modify 2015-03-25
	//bool quit_mesia = false;
	bool quit_abacus = false;

	ifstream ifs;
	if(MY_RANK==0)
	{
		ifs.open(fn.c_str());
	 	if (!ifs)
	 	{
			//xiaohui modify 2015-03-25
			//quit_mesia = true;
			quit_abacus = true;
	 	}
		else
		{
			// if the number is not match,
			// quit the program or not.
			bool quit=false;

			string name;
			ifs >> name;

			// check lattice constant, unit is Angstrom
			CHECK_DOUBLE(ifs,ucell.lat0 * BOHR_TO_A,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e11,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e12,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e13,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e21,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e22,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e23,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e31,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e32,quit);
			CHECK_DOUBLE(ifs,ucell.latvec.e33,quit);

			for(int it=0; it<ucell.ntype; it++)
			{
				CHECK_STRING(ifs,ucell.atoms[it].label,quit);
			}

			for(int it=0; it<ucell.ntype; it++)
			{
				CHECK_DOUBLE(ifs,ucell.atoms[it].na,quit);
			}

			string coordinate;
			ifs >> coordinate;

			for(int it=0; it<ucell.ntype; it++)
			{
				for(int ia=0; ia<ucell.atoms[it].na; ia++)
				{
					CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].x,quit);
					CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].y,quit);
					CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].z,quit);
				}
			}

			CHECK_INT(ifs, NSPIN);
			if(NSPIN == 1)
			{
				READ_VALUE(ifs, en.ef);
				ofs_running << " read in fermi energy = " << en.ef << endl;
			}
			else if(NSPIN == 2)
			{
				if(is==0)READ_VALUE(ifs, en.ef_up);
				else if(is==1)READ_VALUE(ifs, en.ef_dw);
			}
			else
			{
				WARNING_QUIT("read_dm","check nspin!");
			}
			CHECK_INT(ifs, NLOCAL);
			CHECK_INT(ifs, NLOCAL);
		}// If file exist, read in data.
	} // Finish reading the first part of density matrix.


#ifndef __MPI
	ofs_running << " Read SPIN = " << is+1 << " density matrix now." << endl;

	if(GAMMA_ONLY_LOCAL)
	{
		for(int i=0; i<NLOCAL; ++i)
		{
			for(int j=0; j<NLOCAL; ++j)
			{
				ifs >> DM[is][i][j];
			}
		}
	}
	else
	{
		WARNING_QUIT("Local_Orbital_Charge::read_dm","The nnrg should not be update");
		CHECK_INT(ifs,LNNR.nnrg);

		for(int i=0; i<LNNR.nnrg; ++i)
		{
			ifs >> DM_R[is][i];
		}
	}
#else

	// distribution of necessary data
	//xiaohui modify 2015-03-25
	//Parallel_Common::bcast_bool(quit_mesia);
	Parallel_Common::bcast_bool(quit_abacus);
	//xiaohui modify 2015-03-25
	//if(quit_mesia)
	if(quit_abacus)
	{
	 	WARNING_QUIT("Local_Orbital_Charge::read_dm","Can not find the density matrix file.");
	}


	if(NSPIN==1)
	{
		Parallel_Common::bcast_double(en.ef);
	}
	else if(NSPIN==2)
	{
		Parallel_Common::bcast_double(en.ef_up);
		Parallel_Common::bcast_double(en.ef_dw);
	}


	if(GAMMA_ONLY_LOCAL)
	{
		//ofs_running << " NLOCAL=" << NLOCAL << endl;
		//ofs_running << " lgd_now=" << lgd_now << endl;
		//ofs_running << " GridT.lgd=" << GridT.lgd << endl;

		double *tmp = new double[NLOCAL];
		for(int i=0; i<NLOCAL; ++i)
		{
			//ofs_running << " i=" << i << endl;
			ZEROS(tmp, NLOCAL);
			if(MY_RANK==0)
			{
				for(int j=0; j<NLOCAL; ++j)
				{
					ifs >> tmp[j];
				}
			}
			Parallel_Common::bcast_double(tmp, NLOCAL);

			const int mu = GridT.trace_lo[i];
			if(mu >= 0)
			{	
				for(int j=0; j<NLOCAL; ++j)
				{
					const int nu = GridT.trace_lo[j];
					if(nu >= 0)
					{
						DM[is][mu][nu] = tmp[j];
					}
				}
			}
		}// i
		delete[] tmp;
	}
	else
	{
		WARNING_QUIT("Local_Orbital_Charge::read_dm","not ready to readin DM_R");
	}
#endif
	if(MY_RANK==0) ifs.close();

	ofs_running << " Finish reading density matrix." << endl;

	timer::tick("Local_Orbital_Charge","read_dm");
	return;
}
