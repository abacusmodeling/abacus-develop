#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "blas_interface.h"

#include "src_external/src_test/test_function.h"
#include "src_external/src_test/src_global/complexmatrix-test.h"

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
	//band_local = new int[1];
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

	wfc_dm_2d.init();		// Peize Lin test 2019-01-16
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
	
	wfc_dm_2d.init();		// Peize Lin test 2019-01-16
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
			wfc_dm_2d.cal_dm();		// Peize Lin test 2019-01-16
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

	for(int ik=0; ik<kv.nks; ik++)
	{
		complex<double> **wfc = LOWF.WFC_K[ik];
		const int ispin = kv.isk[ik];
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
				const double wg_local=wf.wg(ik,ib);
				if(wg_local>0)
				{
					if(nRow==0) ibStart=ib;
					const int iline=nRow*nw1;
					complex<double> phase=exp_R*wg_local;
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

//added by zhengdy-soc, for non-collinear case
inline void cal_DM_ATOM_nc(const Grid_Technique &gt, const complex<double> fac, Record_adj RA,
				   const int ia1, const int iw1_lo, const int nw1, const int gstart, 
				   complex<double> *WFC_PHASE, complex<double> **DM_ATOM)
{
	if(NSPIN !=4 ) WARNING_QUIT("Local_Orbital_Charge","NSPIN not match!");
	const char transa='N', transb='T';	
	const complex<double> alpha=1, beta=1;
	int ispin=0;

	for(int is1=0;is1<2;is1++)
	{
		for(int is2=0;is2<2;is2++)
		{
			for(int ik=0; ik<kv.nks; ik++)
			{
				complex<double> **wfc = LOWF.WFC_K[ik];
				int atom2start=0;

				for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
				{
					complex<double> *DM=&DM_ATOM[ispin][atom2start];
					const int T2 = RA.info[ia1][ia2][3];
					const int I2 = RA.info[ia1][ia2][4];
					Atom* atom2 = &ucell.atoms[T2];
					const int start2 = ucell.itiaiw2iwt(T2,I2,0);
					const int iw2_lo=gt.trace_lo[start2] + NLOCAL/NPOL*is2;
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
								WFC_PHASE[iline+iw1]=phase*conj(wfc[ib][iw1_lo+iw1 + NLOCAL/NPOL*is1]);
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
			ispin++;
		}//is2
	}//is1
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
				if(!NONCOLIN)cal_DM_ATOM(gt, fac, RA, ca, iw1_lo, nw1, gstart, WFC_PHASE, DM_ATOM);
				else cal_DM_ATOM_nc(gt, fac, RA, ca, iw1_lo, nw1, gstart, WFC_PHASE, DM_ATOM);

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
								const double wg_local = wf.wg(is, ib);
								if(wg_local>0)
								{
									// dm = \sum ( wg * c[ib][mu] * c[ib][nu] )
									// dm saved accordint to sub-FFT box.
									alpha[nu_local] += wg_local * (conj(LOWF.WFC_GAMMA_B[is][ib][mu_local]) * LOWF.WFC_GAMMA_B[is][ib][nu_local]).real();
								}
							}
						}
					}
				}
			}
		}
	}
#ifdef __MPI //2015-09-06, xiaohui
	else	// Peize Lin update 2018-07-02
	{	
		for( int is=0; is<NSPIN; ++is )
			for (int i=0; i<lgd_now; i++)
				ZEROS(this->DM[is][i], lgd_now);

		int nprocs,myid;
		MPI_Status status;
		MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
		MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

		vector<int> bands_local(DSIZE);
		for (int id=0; id<DSIZE; id++)
			bands_local[id] = (id<NBANDS%DSIZE) ? NBANDS/DSIZE+1 : NBANDS/DSIZE;
		const int band_local = bands_local[DRANK];
		
		int lastband_in_proc = 0;
		for (int id=0, count_bands=0; id<DSIZE; id++)
		{
			count_bands += bands_local[id];
			if (count_bands >= NBANDS)
			{
				lastband_in_proc = id;
				break;
			}
		}
		
		matrix wg_local(NSPIN,band_local);
		for(int id=0, Total_Bands=0; id <= lastband_in_proc; ++id)
		{
			if(myid == id)
				for(int is=0; is<NSPIN; is++)
					for (int ib=0; ib<bands_local[myid]; ib++)
						wg_local(is,ib) = wf.wg(is,Total_Bands+ib);
			Total_Bands += bands_local[id];
		}

		for( int is=0; is<NSPIN; ++is )
		{
			matrix Z_wg( NLOCAL, band_local );
			if(myid <= lastband_in_proc)
				for(int iw=0; iw<NLOCAL; iw++)
					for(int ib=0; ib<bands_local[myid]; ib++)
						Z_wg(iw,ib) = ParaO.Z_LOC[is][iw*bands_local[myid]+ib] * wg_local(is,ib);

			const int row_col = (NLOCAL%300) ? NLOCAL/300+1 : NLOCAL/300;
				
			matrix Z_row;
			matrix Z_col;
			matrix rho_row_col;

			for(int row_count=0; row_count<row_col; row_count++)
			{
				const int row_remain = ( (row_count+1)*300 <= NLOCAL )
									 ? 300
									 : NLOCAL - row_count*300;
				
				Z_row.create( row_remain, band_local, false );
				for(int i_row=0; i_row<row_remain; i_row++)
				{
					const int row_index = row_count*300 + i_row;
					for(int ib=0; ib<band_local; ib++)
						Z_row(i_row,ib) = Z_wg(row_index,ib);
				}

				for(int col_count=0; col_count<row_col; col_count++)
				{
					const int col_remain = ( (col_count+1)*300 <= NLOCAL )
										 ? 300
										 : NLOCAL - col_count*300;
												
					Z_col.create( col_remain, band_local, false );
					for(int i_col=0; i_col<col_remain; i_col++)
					{
						const int col_index = i_col +col_count*300;
						for(int ib=0; ib<band_local; ib++)
							Z_col(i_col,ib) = ParaO.Z_LOC[is][col_index*band_local+ib] ;
					}
					
					rho_row_col.create( row_remain, col_remain, false );
					
					//for(int i_row=0; i_row<row_remain; i_row++)
					//	for(int i_col=0; i_col<col_remain; i_col++)
					//		for(int ib=0; ib<band_local; ib++)
					//			rho_row_col(i_row,i_col) += Z_row(i_row,ib) * Z_col(i_col,ib);
										
					LapackConnector::gemm(
						'N', 'T', 
						row_remain, col_remain, band_local,
						1, Z_row.c, band_local, Z_col.c, band_local,
						0, rho_row_col.c, col_remain);
					MPI_Barrier(DIAG_HPSEPS_WORLD);
					Parallel_Reduce::reduce_double_all( rho_row_col.c, row_remain*col_remain);

					if(GAMMA_ONLY_LOCAL)
					{
						for(int i_row=0; i_row<row_remain; i_row++)
						{
							const int row_index = row_count*300 + i_row;
							const int row_mu = GridT.trace_lo[row_index];
							if(row_mu<0)	continue;
							for(int i_col=0; i_col<col_remain; i_col++)
							{
								const int col_index = col_count*300 + i_col;
								const int col_nu = GridT.trace_lo[col_index];
								if(col_nu<0)	continue;
								this->DM[is][row_mu][col_nu] = rho_row_col(i_row,i_col);
							}
						}
					}
				}  // end for col_count
			}  // end for row_count
		}  // end for is	
	}  // end if !BFIELD
#endif //2015-09-06, xiaohui
#ifndef __MPI //2015-09-06, xiaohui
	else
	{
		for(int is=0; is<NSPIN; is++)
			for (int i=0; i<lgd_now; i++)
				ZEROS(this->DM[is][i], lgd_now);
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
								const double wg_local = wf.wg(is, ib);
								//ofs_running << " wg_local=" << wg_local << endl;
								if(wg_local>0)
								{
									// dm = \sum ( wg * c[ib][mu] * c[ib][nu] )
									// dm saved accordint to sub-FFT box.
									alpha[nu_local] += wg_local * LOWF.WFC_GAMMA[is][ib][mu_local] 
									* LOWF.WFC_GAMMA[is][ib][nu_local];

								}//wg_local
							}//ib
						}//nu_local
					}//j
				}//mu_local
			}//i
		}//is
	}
#endif //2015-09-06, xiaohui
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
		if(NSPIN==1||NSPIN==4)
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
			if(NSPIN == 1||NSPIN == 4)
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


	if(NSPIN==1||NSPIN==4)
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
