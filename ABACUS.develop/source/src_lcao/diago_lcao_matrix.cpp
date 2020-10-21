#include "diago_lcao_matrix.h"
#include "../src_pw/algorithms.h"
#include "../src_pw/global.h"
#include "../src_external/src_pdiag/pdiag_double.h"
#include "../src_lcao/hs_matrix.h"
//xiaohui modified 2013-03-23
//#include "../src_develop/src_siao/selinv.h"

Diago_LCAO_Matrix::Diago_LCAO_Matrix(){}
Diago_LCAO_Matrix::~Diago_LCAO_Matrix(){}

void Diago_LCAO_Matrix::solve_complex_matrix(const int &ik, complex<double>** wfc, ComplexMatrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","solve_complex_matrix");
	time_t time_start = time(NULL);
//	ofs_running << " Start Time : " << ctime(&time_start);
//

	//if(DIAGO_TYPE=="lapack") xiaohui modify 2013-09-02
	if(KS_SOLVER=="lapack") //xiaohui add 2013-09-02
	{
		this->using_LAPACK_complex(ik, wfc);
	}
	else
	{
#ifdef __MPI
		this->using_HPSEPS_complex(ik, wfc, wfc_2d);
#else
		WARNING_QUIT("Diago_LCAO_Matrix::solve_complex_matrix","only lapack is available!");
#endif
	}

	time_t time_end = time(NULL);
//	ofs_running << " End   Time : " << ctime(&time_end);
//	ofs_running << " FINAL Time : " << difftime(time_end, time_start) << " (SEC)" << endl;

	OUT_TIME("diago(complex)", time_start, time_end);

	return;
}

void Diago_LCAO_Matrix::solve_double_matrix(const int &ik, double** wfc, matrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","solve_double_matrix");
	timer::tick("Diago_LCAO_Matrix","solve_double_matrix",'F');
	time_t time_start = time(NULL);
//	ofs_running << " Start Time : " << ctime(&time_start);

	// mohan test the combination of hpseps and cg.


	//DIAGO_TYPE="hpseps";
	/*
	static int count = 1;
	if(count < 2)
	{
		cout << " first step is hpseps! " << endl;
		DIAGO_TYPE="hpseps";
	}
	else
	{
		DIAGO_TYPE="cg";
	}
	++count;
	*/
	
	//if(DIAGO_TYPE=="lapack") xiaohui modify 2013-09-02
	if(KS_SOLVER=="lapack") //xiaohui add 2013-09-02
	{
		this->using_LAPACK(ik, wfc);
	}
/*
#ifdef __SELINV
	//xiaohui modified 2013-03-23
	else if(DIAGO_TYPE=="selinv") //mohan add 2011-09
	{
		Selinv SV;
		// job=1: density matrix for charge density.
		// job=2: density matrix for force.
		// job=3: energy density matrix for force.
		int job = 1;
		SV.using_SELINV(ik, job, LM.Hloc, LM.Sloc);
	}
#endif
	else if(DIAGO_TYPE=="cg")
	{
		this->using_CG(ik, wfc);
	}
*/
#ifdef __MPI
	//else if(DIAGO_TYPE=="hpseps") xiaohui modify 2013-09-02
	else if(KS_SOLVER=="hpseps" || KS_SOLVER=="genelpa") //yshen add 7/15/2016
	{
		this->using_HPSEPS_double(ik, wfc, wfc_2d);
	}
#endif
	else
	{
		//cout << " Diago_LCAO_Matrix, diago_type = " << DIAGO_TYPE << endl; xiaohui modify 2013-09-02
		cout << " Diago_LCAO_Matrix, diago_type = " << KS_SOLVER << endl; //xiaohui add 2013-09-02
		//WARNING_QUIT("Diago_LCAO_Matrix::init","Check DIAGO_TYPE."); xiaohui modify 2013-09-02
		WARNING_QUIT("Diago_LCAO_Matrix::init","Check KS_SOLVER."); //xiaohui add 2013-09-02
	}
	time_t time_end = time(NULL);
//	ofs_running << " End   Time : " << ctime(&time_end);
//	ofs_running << " FINAL Time : " << difftime(time_end, time_start) << " (SEC)" << endl;

	OUT_TIME("diago(double)",time_start, time_end);

	/*
	ofs_running << "\n Output c:" << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
//			ofs_running << " " << setw(8) << i << setw(8) << j << setw(25) << c[i][j] << endl;
			if( abs(c[i][j]) > 1.0e-5 )
			{ 
				ofs_running << setw(20) << c[i][j];		
			}
			else
			{
				ofs_running << setw(20) << "0";
			}	
		}
		ofs_running << endl;
	}
	*/

	timer::tick("Diago_LCAO_Matrix","solve_double_matrix",'F');
	return;
}

#ifdef __MPI
void Diago_LCAO_Matrix::using_HPSEPS_double(const int &ik, double**wfc, matrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","using_HPSEPS_double");

	// save H and S matrix to disk.
	bool bit = false;
	HS_Matrix::saving_HS(LM.Hloc, LM.Sloc, bit, ParaO.out_hs);
	ofs_running << setprecision(6);

	// Distribution of matrix for 
	// prallel eigensolver.
	ParaO.diago_double_begin(ik, wfc, wfc_2d, LM.Hloc, LM.Sloc, wf.ekb[ik]);

	/*
	string fh = "/home/mohan/3_my_program/1_DAPE/data/data-H32";
	string fs = "/home/mohan/3_my_program/1_DAPE/data/data-S32";
	double* eigen = new double[32];
	double* eigvr = new double[32*32];
	PD.readin(fh, fs, 32, eigen, eigvr);
	*/

	return;
}


void Diago_LCAO_Matrix::using_HPSEPS_complex(const int &ik, complex<double>** wfc, ComplexMatrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","using_HPSEPS_complex");

	//ParaO.out_hs=1;//zhengdy-soc-test
	bool bit = false; //LiuXh, 2017-03-21
	//if set bit = true, there would be error in soc-multi-core calculation, noted by zhengdy-soc
	HS_Matrix::saving_HS_complex(LM.Hloc2, LM.Sloc2, bit, ParaO.out_hs); //LiuXh, 2017-03-21
	ofs_running << setprecision(6); //LiuXh, 2017-03-21

	ParaO.diago_complex_begin(ik, wfc, wfc_2d, LM.Hloc2, LM.Sloc2, wf.ekb[ik]);
	//added by zhengdy-soc, rearrange the WFC_K from [up,down,up,down...] to [up,up...down,down...], 
	if(NONCOLIN)
	{
		int row = GridT.lgd;
		vector<complex<double>> tmp(row);
		for(int ib=0; ib<NBANDS; ib++)
		{
			for(int iw=0; iw<row / NPOL; iw++)
			{
				tmp[iw] = LOWF.WFC_K[ik][ib][iw * NPOL];
				tmp[iw + row / NPOL] = LOWF.WFC_K[ik][ib][iw * NPOL + 1];
			}
			for(int iw=0; iw<row; iw++)
			{
				LOWF.WFC_K[ik][ib][iw] = tmp[iw];
			}
		}
	}

	return;
}
#endif

void Diago_LCAO_Matrix::using_LAPACK_complex(const int &ik, complex<double> **wfc)const
{
	TITLE("Diago_LCAO_Matrix","using_LAPACK_complex");

	assert(NPROC = 1);

	ComplexMatrix Htmp(NLOCAL,NLOCAL);
	ComplexMatrix Stmp(NLOCAL,NLOCAL);
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
//			cout<<"H_and_S: "<<i<<" "<<j<<" "<<LM.Hloc2[i*NLOCAL+j]<<" "<<LM.Sloc2[i*NLOCAL+j]<<endl;
			Htmp(i,j) = LM.Hloc2[i*NLOCAL+j];
			Stmp(i,j) = LM.Sloc2[i*NLOCAL+j];
		}
	}

//----------------------------
// keep this for tests
//    out.printcm_norm("Lapack_H", Htmp, 1.0e-5);
//    out.printcm_norm("Lapack_S", Stmp, 1.0e-5);
//----------------------------

	double* en = new double[NLOCAL];
	ZEROS(en, NLOCAL);

	ComplexMatrix hvec(NLOCAL, NBANDS);
	hm.cdiaghg(NLOCAL, NBANDS, Htmp, Stmp, NLOCAL, en, hvec);

	if(!NONCOLIN)
	for(int ib=0; ib<NBANDS; ib++)
	{
		for(int iw=0; iw<NLOCAL; iw++)
		{
			LOWF.WFC_K[ik][ib][iw] = hvec(iw,ib);
		}
	}
	else
	for(int ib=0; ib<NBANDS; ib++)
	{
		for(int iw=0; iw<NLOCAL / NPOL; iw++)
		{
			LOWF.WFC_K[ik][ib][iw] = hvec(iw * NPOL, ib);
			LOWF.WFC_K[ik][ib][iw + NLOCAL / NPOL] = hvec(iw * NPOL + 1, ib);
		}
	}

	//cout << "\n Energy for k=" << ik << endl; 
	for(int ib=0; ib<NBANDS; ib++)
	{
		wf.ekb[ik][ib] = en[ib];
		//cout << setw(10) << ib << setw(15) << wf.ekb[ik][ib] * Ry_to_eV << endl;
	}

	return;
}

void Diago_LCAO_Matrix::using_LAPACK(const int &ik, double** wfc)const
{
	TITLE("Diago_LCAO_Matrix","using_LAPACK");
	assert(NLOCAL>0);

	// save H and S matrix to disk.
//	bool bit = false;
	bool bit = true;//zhengdy-soc
	HS_Matrix::saving_HS(LM.Hloc, LM.Sloc, bit, ParaO.out_hs);

	matrix Htmp(NLOCAL,NLOCAL);
	matrix Stmp(NLOCAL,NLOCAL);
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			Htmp(i,j) = LM.Hloc[i*NLOCAL+j];
			Stmp(i,j) = LM.Sloc[i*NLOCAL+j];
		}
	}

	// @@@@@@@
	// test
	// @@@@@@@
//    out.printrm("Lapack_H", Htmp);
//    out.printrm("Lapack_S", Stmp);
		
	int itype=1;
	int lwork=3*NLOCAL-1;// tmp
	double* w = new double[NLOCAL];
	double* work = new double[lwork];
	ZEROS(w, NLOCAL);
	ZEROS(work, lwork);
	int info;

	clock_t clock_start, clock_end;
	clock_start = std::clock();
	LapackConnector::dsygv(itype,'V','U',NLOCAL,Htmp,NLOCAL,Stmp,NLOCAL,w,work,lwork,&info);
	clock_end = std::clock();
	double duration = (double)(clock_end - clock_start) / CLOCKS_PER_SEC;

	ofs_running << setiosflags(ios::fixed) << setprecision(20);
//	ofs_running << " clock_start = " << clock_start << endl;
//	ofs_running << " clock_end = " << clock_end << endl;
	ofs_running << " Time using dsygv in LAPACK (seconds) is " << duration << endl;
	ofs_running << resetiosflags(ios::fixed) << setprecision(10);

	for(int i=0; i<NBANDS; i++)
	{
		// eigenvalues
		wf.ekb[ik][i] = w[i]; 
		for(int j=0; j<NLOCAL; j++)
		{
			wfc[i][j] = Htmp(j,i);
		}
	}
	

	// @@@@@@@
	// test
	// @@@@@@@
	/*
	cout << "\n Lapack, wfc after diago:" << endl;
	for(int i=0; i<NBANDS; i++)
	{
		cout << " Eigenvalue from LAPACK : " << setw(5) << setw(12) << wf.ekb[ik][i] << endl;
		cout << " Eigenfunctions" << endl;
		for(int j=0; j<NLOCAL; j++)
		{
			cout << setw(12) << wfc[i][j];
		}
		cout << endl;
	}
	//exit(0);
	*/

	delete[] w;
	delete[] work;
	return;
}


#include "local_orbital_elec.h"
//#include "../src_develop/src_cg/cg_lcao.h"
//#include "../src_develop/src_cg/cg_precon.h"
/*
void Diago_LCAO_Matrix::using_CG(const int &ik, double **c)const
{
	TITLE("Diago_LCAO_Matrix","using_CG");
	timer::tick("Diago_LCAO_Matrix","using_CG",'I');
	int notconv = 0;
	double cg_iter=0.0;

	//  out.printcm_norm("Cij_0", this->Cij, 1.0e-5);
	//  ZEROS(en, NBANDS);

//	cout << " ParaA.nlocal=" << ParaA.nlocal << endl;
//	cout << " ParaO.nrow=" << ParaO.nrow << endl;
	
	double *precon = new double[ParaA.nlocal];
	assert(ParaA.nlocal == ParaO.nrow);
	CG_precon::select_precondition(DIAGO_CG_PREC, precon);


	bool reorder = true;
	//bool reorder = false;
	if(ATOM_DISTRIBUTION)
	{
		// 'c' is the wave functions after parallized.
		// because c is adapted to grid integration,
		// and used to generate density matrix,
		// however, ccg is here is divided by atoms and adapted
		// to the CG algorithm.
		double** ccg = new double*[NBANDS];
		for(int ib=0; ib<NBANDS; ++ib)
		{
			ccg[ib] = new double[ParaA.nlocal];
			ZEROS(ccg[ib], ParaA.nlocal);
		}

		// buffer contains the whole wave functions,
		// and is put into ccg.
		double *buffer = new double[NLOCAL];
		int *countb = new int[NLOCAL];
		for(int ib=0; ib<NBANDS; ++ib)
		{
			// put the wave functions in buffer.
			ZEROS(buffer, NLOCAL);		
			for(int iw=0; iw<NLOCAL; ++iw)
			{
				const int nu = GridT.trace_lo[iw];
				if(nu>=0)
				{ 
					buffer[iw] = c[ib][nu];
				}
			}
			Parallel_Reduce::reduce_double_all(buffer,NLOCAL);

			// collect the trace_lo;
			ZEROS(countb, NLOCAL);
			for(int iw=0; iw<NLOCAL; ++iw)
			{
				const int nu = GridT.trace_lo[iw];
				if(nu>=0)
				{
					countb[iw]=1;
				}
			}
			// collect countb 
			Parallel_Reduce::reduce_int_all(countb, NLOCAL);
			// and then get the correct final wave functions
			for(int iw=0; iw<NLOCAL; ++iw)
			{
				assert(countb[iw]>0);
				buffer[iw] /= (double)countb[iw];
			}

			// transfer wave functions from 'c' to ccg
			for(int iw=0; iw<NLOCAL; ++iw)
			{
				const int mu = ParaO.trace_loc_row[iw];	
				if( mu>=0 )
				{
					ccg[ib][mu] = buffer[iw];
				}
			}
		}
		delete[] countb;
		 
		CG_LCAO CGL;
		CGL.diag( ccg, wf.ekb[ik], ParaA.nlocal, NBANDS, precon, ETHR, DIAGO_CG_MAXITER, reorder, notconv, cg_iter);

		for(int ib=0; ib<NBANDS; ++ib)
		{
			// put the local wave function into a big array.
			ZEROS(buffer,NLOCAL);
			for(int iw=0; iw<NLOCAL; ++iw)
			{
				const int mu = ParaO.trace_loc_row[iw];
				if( mu >= 0 )
				{
					buffer[iw] = ccg[ib][mu];	
				}
			}
			// reduce the wave functions. 
			Parallel_Reduce::reduce_double_all(buffer, NLOCAL);	

			// put the wave functions back to LOWF.WFC_GAMMA 
			for(int iw=0; iw<NLOCAL; ++iw)
			{
				const int mu = GridT.trace_lo[iw];
				if(mu >= 0)
				{
					c[ib][mu] = buffer[iw];
				}


				const int mu2 = LOWF.trace_aug[iw];
				if(GAMMA_ONLY_LOCAL)
				{
					if(mu2 >= 0)
					{
						LOWF.WFC_GAMMA_aug[CURRENT_SPIN][ib][mu2] = buffer[iw];
					}
				}
				else
				{
					WARNING_QUIT("using_CG","how to update augmented wave functions.");
				}
			}
		}	
		
		delete[] buffer;

		for(int ib=0; ib<NBANDS; ++ib)
		{
			delete[] ccg[ib];
		}
		delete[] ccg;
	}
	else
	{
		Diago_CG_Real cgr;
		cgr.diag( c, wf.ekb[ik], NLOCAL, NBANDS, precon, ETHR, DIAGO_CG_MAXITER, reorder, notconv, cg_iter);
	}

	delete[] precon;

	//OUT(ofs_running,"Sparse H occupied",hm.hon.UHM.Sparse_H.rate());
	OUT(ofs_running,"cg_iter",cg_iter);
	Local_Orbital_Elec::avg_iter = cg_iter;
	
	//  out.printcm_norm("Cij", this->Cij, 1.0e-5);
	//  out.printr1_d("en",en,NBANDS);

	timer::tick("Diago_LCAO_Matrix","using_CG",'i');
	return;	
}
*/
