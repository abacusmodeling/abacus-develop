#include "force_lcao_gamma.h"
#include "../src_pw/global.h"

Force_LCAO_gamma::Force_LCAO_gamma ()
{
}

Force_LCAO_gamma::~Force_LCAO_gamma ()
{
}

// be called in force_lo.cpp
void Force_LCAO_gamma::ftable_gamma (void)
{
    TITLE("Force_LCAO_gamma", "ftable");
	timer::tick("Force_LCAO_gamma","ftable_gamma",'F');
	
	// allocate DSloc_x, DSloc_y, DSloc_z
	// allocate DHloc_fixed_x, DHloc_fixed_y, DHloc_fixed_z
	this->allocate_gamma();

	// calculate the 'energy density matrix' here.
	this->cal_foverlap();

	// calculate the 'density matrix' here.
	double** dm2d = new double*[NSPIN];
	for(int is=0; is<NSPIN; ++is)
	{
		dm2d[is] = new double[ParaO.nloc];
		ZEROS(dm2d[is], ParaO.nloc);
	}
	Memory::record ("Force_LCAO_gamma", "dm2d", ParaO.nloc*NSPIN, "double");	


	bool with_energy = false;

	for(int is=0; is<NSPIN; ++is)
	{
		this->set_EDM_gamma(dm2d[is], with_energy, is);
	}

	this->cal_ftvnl_dphi(dm2d);
	this->cal_fvnl_dbeta(dm2d);

	// calculate < dphi | V | phi > on real space grid.
	this->cal_fvl_dphi(dm2d);
               

	for(int is=0; is<NSPIN; ++is)
	{
		delete[] dm2d[is];
	}
	delete[] dm2d;

    for (int iat=0; iat<ucell.nat; iat++)
    {
        Parallel_Reduce::reduce_double_pool( this->foverlap[iat], 3);
        Parallel_Reduce::reduce_double_pool( this->ftvnl_dphi[iat], 3);
        Parallel_Reduce::reduce_double_pool( this->fvnl_dbeta[iat], 3);
        Parallel_Reduce::reduce_double_pool( this->fvl_dphi[iat], 3);
    }
	if(STRESS)for (int ipol=0; ipol<3; ipol++)
	{
		Parallel_Reduce::reduce_double_pool( this->soverlap[ipol], 3);
		Parallel_Reduce::reduce_double_pool( this->stvnl_dphi[ipol], 3);
		Parallel_Reduce::reduce_double_pool( this->svnl_dbeta[ipol], 3);
		Parallel_Reduce::reduce_double_pool( this->svl_dphi[ipol], 3);
	}

	// delete DSloc_x, DSloc_y, DSloc_z
	// delete DHloc_fixed_x, DHloc_fixed_y, DHloc_fixed_z
	this->finish_ftable_gamma();

	timer::tick("Force_LCAO_gamma","ftable_gamma",'F');
    return;
}

void Force_LCAO_gamma::allocate_gamma(void)
{
	TITLE("Force_LCAO_gamma","allocate_gamma");
	timer::tick("Force_LCAO_gamma","allocate_gamma");

	// need to calculate the derivative in build_ST_new
	bool cal_deri = true;

    //calculate dS in LCAO
    //liaochen add on 2010/7/12
	//save the results in dense matrix by now.
	//ParaO.nloc: number of H elements in this proc.
    LM.DSloc_x = new double [ParaO.nloc];
    LM.DSloc_y = new double [ParaO.nloc];
    LM.DSloc_z = new double [ParaO.nloc];
    ZEROS(LM.DSloc_x, ParaO.nloc);
    ZEROS(LM.DSloc_y, ParaO.nloc);
    ZEROS(LM.DSloc_z, ParaO.nloc);
	//allocate stress part in gamma_only-line, added by zhengdy-stress
if(STRESS){
	LM.DSloc_11 = new double [ParaO.nloc];
	LM.DSloc_12 = new double [ParaO.nloc];
	LM.DSloc_13 = new double [ParaO.nloc];
	LM.DSloc_22 = new double [ParaO.nloc];
	LM.DSloc_23 = new double [ParaO.nloc];
	LM.DSloc_33 = new double [ParaO.nloc];
	ZEROS(LM.DSloc_11, ParaO.nloc);
	ZEROS(LM.DSloc_12, ParaO.nloc);
	ZEROS(LM.DSloc_13, ParaO.nloc);
	ZEROS(LM.DSloc_22, ParaO.nloc);
	ZEROS(LM.DSloc_23, ParaO.nloc);
	ZEROS(LM.DSloc_33, ParaO.nloc);
	LM.DHloc_fixed_11 = new double [ParaO.nloc];
	LM.DHloc_fixed_12 = new double [ParaO.nloc];
	LM.DHloc_fixed_13 = new double [ParaO.nloc];
	LM.DHloc_fixed_22 = new double [ParaO.nloc];
	LM.DHloc_fixed_23 = new double [ParaO.nloc];
	LM.DHloc_fixed_33 = new double [ParaO.nloc];
	ZEROS (LM.DHloc_fixed_11, ParaO.nloc);
	ZEROS (LM.DHloc_fixed_12, ParaO.nloc);
	ZEROS (LM.DHloc_fixed_13, ParaO.nloc);
	ZEROS (LM.DHloc_fixed_22, ParaO.nloc);
	ZEROS (LM.DHloc_fixed_23, ParaO.nloc);
	ZEROS (LM.DHloc_fixed_33, ParaO.nloc);
}
    //calculate dS in LCAO basis
    // tips: build_ST_new --> ParaO.set_force 
	UHM.UOM.build_ST_new ('S', cal_deri);

	Memory::record("force_lo", "dS", ParaO.nloc*3, "double");

    //calculate dT in LCAP
    //allocation dt
    //liaochen add on 2010/7/12
    LM.DHloc_fixed_x = new double [ParaO.nloc];
    LM.DHloc_fixed_y = new double [ParaO.nloc];
    LM.DHloc_fixed_z = new double [ParaO.nloc];
    ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_z, ParaO.nloc);
    
	//calculate dT
    //calculate T + VNL(P1) in LCAO basis
    UHM.UOM.build_ST_new ('T', cal_deri);
	//test_gamma(LM.DHloc_fixed_x, "dHloc_fixed_x T part");
    
	//UHM.UOM.build_Nonlocal_beta (cal_deri);
    UHM.UOM.build_Nonlocal_mu (cal_deri);
	//test_gamma(LM.DHloc_fixed_x, "dHloc_fixed_x Vnl part");

	Memory::record("force_lo", "dTVNL", ParaO.nloc*3, "double");

	timer::tick("Force_LCAO_gamma","allocate_gamma");
	return;
}

void Force_LCAO_gamma::finish_ftable_gamma(void)
{
    delete [] LM.DSloc_x;
    delete [] LM.DSloc_y;
    delete [] LM.DSloc_z;
    delete [] LM.DHloc_fixed_x;
    delete [] LM.DHloc_fixed_y;
    delete [] LM.DHloc_fixed_z;
	if(STRESS)//added by zhengdy-stress
	{
		delete [] LM.DSloc_11;
		delete [] LM.DSloc_12;
		delete [] LM.DSloc_13;
		delete [] LM.DHloc_fixed_11;
		delete [] LM.DHloc_fixed_12;
		delete [] LM.DHloc_fixed_13;
		delete [] LM.DSloc_22;
		delete [] LM.DSloc_23;
		delete [] LM.DSloc_33;
		delete [] LM.DHloc_fixed_22;
		delete [] LM.DHloc_fixed_23;
		delete [] LM.DHloc_fixed_33;
	}
	return;
}

double Force_LCAO_gamma::set_EDM_element(
	const int &ii, const int &jj, 
	const bool with_energy,
	double*** coef1, double*** coef2, const int &is)
{
	double ene = 0.0;
	for (int ib = 0; ib < NBANDS; ib++)
	{
		const double w1 = wf.wg(is,ib);
		if(w1>0.0)
		{
			if(with_energy)
			{
				ene += w1 * wf.ekb[is][ib] * coef1[is][ib][ii] * coef2[is][ib][jj];
			}
			else
			{
				ene += w1 * coef1[is][ib][ii] * coef2[is][ib][jj];
			}
		}
	}
	return ene;
}

//#include "../src_develop/src_siao/selinv.h"
void Force_LCAO_gamma::set_EDM_gamma(double* dm, bool with_energy, const int &ispin)//ispin is the spin index. 
{
	TITLE("Force_LCAO_gamma","set_EDM_gamma");
	timer::tick("Force_LCAO_gamma","set_EDM");

#ifdef __SELINV
//xiaohui modified 2013-03-23, adding "/*"
/*    if(DIAGO_TYPE=="selinv")
	{
//		cout << " set_EDM_gamma()" << endl;
//		cout << " fill dm with density matrix or energy density matrix here." << endl;

		//--------------------------------------------------
		// job=1: density matrix for charge density (grid).
		// job=2: density matrix for force (dm2d).
		// job=3: energy density matrix for force (dm2d).
		//--------------------------------------------------
		const int ik = 0;
		int job = 0;
		if(with_energy) 
		{
			job=3;
		}
		else 
		{
			job=2;
		}
		Selinv::dm_ptr = dm;
		Selinv::using_SELINV(ik, job, LM.Hloc, LM.Sloc);	

        return;
	}
*/
#endif
	//xiaohui add 2015-03-24, like density matrix calculation
#ifdef __MPI //2015-09-06, xiaohui
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
			for (int ib=0; ib< Nbands; ib++)
			{
				if(with_energy)
				{
					w1[ispin][ib] = wf.wg(ispin, ib + Total_Bands) * wf.ekb[ispin][ib + Total_Bands];
				}
				else
				{
					w1[ispin][ib] = wf.wg(ispin, ib + Total_Bands);
				}
			}
		}
		Total_Bands += Nbands;
	}

	double** Z_wg;
	Z_wg = new double*[NSPIN];
	for(int is=0; is<NSPIN; is++)
	{
		Z_wg[is] = new double[local_band * NLOCAL];
		ZEROS(Z_wg[is], local_band * NLOCAL);
	}

	time_t time_Z_wg_start = time(NULL);
	double Time_Z_Wg_Start = (double)clock();

	if(myid <= lastband_in_proc)
	{
		for(int i=0; i<local_bands[myid]; i++)
		{
			for(int n=0; n<NLOCAL; n++)
			{
				Z_wg[ispin][i + n*local_bands[myid]] = ParaO.Z_LOC[ispin][i + n*local_bands[myid]]* w1[ispin][i];
			}
		}
	}

	for(int is=0; is<NSPIN; is++)
	{
		delete[] w1[is];
	}
	delete[] w1;
	delete[] local_bands;

	time_t time_Z_wg_end = time(NULL);
	double Time_Z_Wg_End =(double)clock();
	double time_Z_wg = difftime(time_Z_wg_start,time_Z_wg_end);

	double coll_time = 0.0;
	double Time_Cal_Start = (double)clock();

	int row_col = NLOCAL/300;
	if(row_col >= 1)
	{
		double** ZW_300;
		ZW_300 = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			ZW_300[is] = new double[300 * local_band];
			ZEROS(ZW_300[is], 300 * local_band);
		}

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
				for(int i_row=0; i_row<300; i_row++)
				{
					row_index = i_row +(row_count-1)*300;
					for(int ib=0; ib<local_band; ib++)
					{
						ZW_300[ispin][ib + i_row*local_band] = Z_wg[ispin][ib + row_index*local_band] ;
					}
				}

				for(int col_count=1; col_count<=row_col; col_count++)
				{
					col_global = col_count*300;
					if(col_global <= NLOCAL)
					{
						double **rho_300;
						rho_300 = new double*[NSPIN];
						for(int is=0; is<NSPIN; is++)
						{
							rho_300[is] = new double[300*300];
							ZEROS(rho_300[is], 300*300);
						}

						for(int i_col=0; i_col<300; i_col++)
						{
							col_index = i_col +(col_count-1)*300;
							for(int ib=0; ib<local_band; ib++)
							{
								Z_300[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index*local_band] ;
							}
						}

					 	double alpha=1;
					 	double beta=0;
					 	char transa = 'T';
					 	char transb = 'N';
					 	int M = 300;
					 	int N = 300;
					 	int K = local_band;
					 	int lda = local_band;
					 	int ldb = local_band;
					 	int ldc = 300;
						dgemm_(&transa, &transb, &M, &N, &K,
							&alpha, &ZW_300[ispin][0], &lda, &Z_300[ispin][0], &ldb,
			 				&beta, &rho_300[ispin][0], &ldc);

						double Time_Coll_Start = (double)clock();
						MPI_Barrier(DIAG_HPSEPS_WORLD);
						Parallel_Reduce::reduce_double_all( rho_300[ispin], 300*300);
						double Time_Coll_End = (double)clock();

						coll_time += (Time_Coll_End-Time_Coll_Start);
						if(GAMMA_ONLY_LOCAL)
        					{
							for(int i_row=0; i_row<300; i_row++)
							{
								row_index = i_row +(row_count-1)*300;
								int row_mu = ParaO.trace_loc_row[row_index];
								for(int i_col=0; i_col<300; i_col++)
								{
									col_index = i_col +(col_count-1)*300;
									int col_nu = ParaO.trace_loc_col[col_index];
									if(row_mu >= 0 && col_nu >= 0)
									{
										int index = row_mu * ParaO.ncol + col_nu;
										//dm[index] = rho_300[ispin][i_col + i_row*300];
										dm[index] = rho_300[ispin][i_row + i_col*300];
										}
								}
							}
						}
						for(int is=0; is<NSPIN; is++)
						{
							delete[] rho_300[is];
						}
						delete[] rho_300;
					}
					else
					{
						int col_remain = 0;
						int col_index2 = 0;
						col_remain = NLOCAL - (col_count - 1)*300;
						double** Z_remain;
						Z_remain = new double*[NSPIN];
						for(int is=0; is<NSPIN; is++)
						{
							Z_remain[is] = new double[col_remain * local_band];
							ZEROS(Z_remain[is], col_remain * local_band);
						}

						double** rho_300_remain;
						rho_300_remain = new double*[NSPIN];
						for(int is=0; is<NSPIN; is++)
						{
							rho_300_remain[is] = new double[300*col_remain];
							ZEROS(rho_300_remain[is], 300*col_remain);
						}
						for(int i_col=0; i_col<col_remain; i_col++)
						{
							col_index2 = i_col +(col_count-1)*300;
							for(int ib=0; ib<local_band; ib++)
							{
								Z_remain[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index2*local_band] ;
							}
						}

					 	double alpha=1;
					 	double beta=0;
					 	char transa = 'T';
					 	char transb = 'N';
					 	int M = 300;
					 	int N = col_remain;
					 	int K = local_band;
					 	int lda = local_band;
					 	int ldb = local_band;
					 	int ldc = 300;
						dgemm_(&transa, &transb, &M, &N, &K,
				  			&alpha, &ZW_300[ispin][0], &lda, &Z_remain[ispin][0], &ldb,
			 					&beta, &rho_300_remain[ispin][0], &ldc);

						double Time_Coll_Start = (double)clock();
						MPI_Barrier(DIAG_HPSEPS_WORLD);
						Parallel_Reduce::reduce_double_all( rho_300_remain[ispin], 300*col_remain);
						double Time_Coll_End = (double)clock();

						coll_time += (Time_Coll_End-Time_Coll_Start);
						if(GAMMA_ONLY_LOCAL)
        					{
							for(int i_row=0; i_row<300; i_row++)
							{
								row_index = i_row +(row_count-1)*300;
								int row_mu = ParaO.trace_loc_row[row_index];
								for(int i_col=0; i_col<col_remain; i_col++)
								{
									col_index = i_col +(col_count-1)*300;
									int col_nu = ParaO.trace_loc_col[col_index];
									if(row_mu >= 0 && col_nu >= 0)
									{
										int index = row_mu * ParaO.ncol + col_nu;
										//dm[index] = rho_300_remain[ispin][i_col + i_row*col_remain];
										dm[index] = rho_300_remain[ispin][i_row + i_col*300];
									}
								}
							}
						}
						for(int is=0; is<NSPIN; is++)
						{
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
				double** Z_row_remain;
				Z_row_remain = new double*[NSPIN];
				for(int is=0; is<NSPIN; is++)
				{
					Z_row_remain[is] = new double[row_remain * local_band];
					ZEROS(Z_row_remain[is], row_remain * local_band);
				}

				for(int i_row=0; i_row<row_remain; i_row++)
				{
					row_index2 = i_row +(row_count-1)*300;
					for(int ib=0; ib<local_band; ib++)
					{
						Z_row_remain[ispin][ib + i_row*local_band] = Z_wg[ispin][ib + row_index2*local_band] ;
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
						}

						for(int i_col=0; i_col<300; i_col++)
						{
							col_index = i_col +(col_count-1)*300;
							for(int ib=0; ib<local_band; ib++)
							{
								Z_300[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index*local_band] ;
							}
						}

					 	double alpha=1;
					 	double beta=0;
					 	char transa = 'T';
					 	char transb = 'N';
					 	int M = row_remain;
					 	int N = 300;
					 	int K = local_band;
					 	int lda = local_band;
					 	int ldb = local_band;
					 	int ldc = row_remain;
						dgemm_(&transa, &transb, &M, &N, &K,
				  			&alpha, &Z_row_remain[ispin][0], &lda, &Z_300[ispin][0], &ldb,
			 				&beta, &rho_remain_300[ispin][0], &ldc);

						double Time_Coll_Start = (double)clock();
						MPI_Barrier(DIAG_HPSEPS_WORLD);
						Parallel_Reduce::reduce_double_all( rho_remain_300[ispin], row_remain*300);
						double Time_Coll_End = (double)clock();

						coll_time += (Time_Coll_End-Time_Coll_Start);
						if(GAMMA_ONLY_LOCAL)
        					{
							for(int i_row=0; i_row<row_remain; i_row++)
							{
								row_index = i_row +(row_count-1)*300;
								int row_mu = ParaO.trace_loc_row[row_index];
								for(int i_col=0; i_col<300; i_col++)
								{
									col_index = i_col +(col_count-1)*300;
									int col_nu = ParaO.trace_loc_col[col_index];
									if(row_mu >= 0 && col_nu >= 0)
									{
										int index = row_mu * ParaO.ncol + col_nu;
										dm[index] = rho_remain_300[ispin][i_row + i_col*row_remain];
									}
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
						double** Z_col_remain;
						Z_col_remain = new double*[NSPIN];
						for(int is=0; is<NSPIN; is++)
						{
							Z_col_remain[is] = new double[col_remain * local_band];
							ZEROS(Z_col_remain[is], col_remain * local_band);
						}
						double** rho_remain_remain;
						rho_remain_remain = new double*[NSPIN];
						for(int is=0; is<NSPIN; is++)
						{
							rho_remain_remain[is] = new double[row_remain * col_remain];
							ZEROS(rho_remain_remain[is], row_remain * col_remain);
						}
						for(int i_col=0; i_col<col_remain; i_col++)
						{
							col_index2 = i_col +(col_count-1)*300;
							for(int ib=0; ib<local_band; ib++)
							{
								Z_col_remain[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index2*local_band] ;
							}
						}

					 	double alpha=1;
					 	double beta=0;
					 	char transa = 'T';
					 	char transb = 'N';
					 	int M = row_remain;
					 	int N = col_remain;
					 	int K = local_band;
					 	int lda = local_band;
					 	int ldb = local_band;
					 	int ldc = row_remain;
						dgemm_(&transa, &transb, &M, &N, &K,
				  			&alpha, &Z_row_remain[ispin][0], &lda, &Z_col_remain[ispin][0], &ldb,
			 				&beta, &rho_remain_remain[ispin][0], &ldc);

						double Time_Coll_Start = (double)clock();
						MPI_Barrier(DIAG_HPSEPS_WORLD);
						Parallel_Reduce::reduce_double_all( rho_remain_remain[ispin], row_remain*col_remain);
						double Time_Coll_End = (double)clock();

						coll_time += (Time_Coll_End-Time_Coll_Start);
						if(GAMMA_ONLY_LOCAL)
        					{
							for(int i_row=0; i_row<row_remain; i_row++)
							{
								row_index = i_row +(row_count-1)*300;
								int row_mu = ParaO.trace_loc_row[row_index];
								for(int i_col=0; i_col<col_remain; i_col++)
								{
									col_index = i_col +(col_count-1)*300;
									int col_nu = ParaO.trace_loc_col[col_index];
									if(row_mu >= 0 && col_nu >= 0)
									{
										int index = row_mu * ParaO.ncol + col_nu;
										dm[index] = rho_remain_remain[ispin][i_row + i_col*row_remain];
									}
								}
							}
						}
						for(int is=0; is<NSPIN; is++)
						{
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
		for(int i_row=0; i_row<NLOCAL; i_row++)
		{
			for(int i_col=0; i_col<NLOCAL; i_col++)
			{
				for(int ib=0; ib<local_band; ib++)
				{
					rho_NLOCAL_NLOCAL[ispin][i_col + i_row*NLOCAL] += Z_wg[ispin][ib + i_row*local_band]  * ParaO.Z_LOC[ispin][ib + i_col*local_band];
				}
			}
		}
		double Time_Coll_Start = (double)clock();
		MPI_Barrier(DIAG_HPSEPS_WORLD);
		Parallel_Reduce::reduce_double_all( rho_NLOCAL_NLOCAL[ispin], NLOCAL*NLOCAL);
		double Time_Coll_End = (double)clock();
		coll_time += (Time_Coll_End-Time_Coll_Start);
		if(GAMMA_ONLY_LOCAL)
        	{
			for(int i_row=0; i_row<NLOCAL; i_row++)
			{
				int row_mu = ParaO.trace_loc_row[i_row];
				for(int i_col=0; i_col<NLOCAL; i_col++)
				{
					int col_nu = ParaO.trace_loc_col[i_col];
					if(row_mu >= 0 && col_nu >= 0)
					{
						int index = row_mu * ParaO.ncol + col_nu;
						dm[index] = rho_NLOCAL_NLOCAL[ispin][i_col + i_row*NLOCAL];
					}
				}
			}
		}
		for(int is=0; is<NSPIN; is++)
		{
			delete[] rho_NLOCAL_NLOCAL[is];
		}
		delete[] rho_NLOCAL_NLOCAL;
	}

	for(int is=0; is<NSPIN; is++)
	{
		delete[] Z_wg[is];
	}
	delete[] Z_wg;
#endif //2015-09-06, xiaohui
#ifndef __MPI //2015-09-06, xiaohui
	// orbital 1
	for(int i=0; i<NLOCAL; i++)
	{
		const int mu = ParaO.trace_loc_row[i];
		if(mu<0) continue;
		const int ii = GridT.trace_lo[i];
		// orbital 2
		for(int j=0; j<NLOCAL; j++)
		{
			const int nu = ParaO.trace_loc_col[j];
			if(nu<0) continue;
			const int jj = GridT.trace_lo[j];
			// energy density matrix
			const int index = mu * ParaO.ncol + nu;
			double ene = 0.0;

			// (1) density matrix can be generated from LOWF.WFC_GAMMA directly.
			if(ii>=0 && jj>=0)
			{
				ene = this->set_EDM_element(ii, jj, with_energy, LOWF.WFC_GAMMA, LOWF.WFC_GAMMA, ispin);
			}

			// (2)
			else if(ii>=0 && jj<0)
			{
				const int a4 = LOWF.trace_aug[j];
				assert(a4>=0);
				ene = this->set_EDM_element(ii, a4, with_energy,  LOWF.WFC_GAMMA,  LOWF.WFC_GAMMA_aug, ispin);
			}
			else if(ii<0 && jj>=0)
			{
				const int a3 = LOWF.trace_aug[i];
				assert(a3>=0);
				// mohan fix serious bug 2011-07-01 (ii, a3) -> (a3, jj) !!!!!!!!!!!!
				ene = this->set_EDM_element(a3, jj, with_energy, LOWF.WFC_GAMMA_aug, LOWF.WFC_GAMMA, ispin);
			}
			else if(ii<0 && jj<0)
			{
				const int a3 = LOWF.trace_aug[i];
				const int a4 = LOWF.trace_aug[j];
				assert(a3>=0);
				assert(a4>=0);
				ene = this->set_EDM_element(a3, a4, with_energy, LOWF.WFC_GAMMA_aug, LOWF.WFC_GAMMA_aug, ispin);
			}

			dm[index] = ene;
			//dm[index] = 1.0;// mohan tmp
		}// end j
	}// end i
#endif //2015-09-06, xiaohui
	timer::tick("Force_LCAO_gamma","set_EDM");
	return;
}


void Force_LCAO_gamma::DerivT_PW(void)
{
    WARNING_QUIT("Force_LCAO_gamma::DerivT_PW","no use for a long time.");
    /*
    //------------------------------------------
    //TEST OUTPUT OF S Deriv Matrix in PW
    //------------------------------------------
    //test under pw basis
    hm.hon.UOM.ForcedT = new matrix[3];
    this->dt = hm.hon.UOM.ForcedT;
    for (int i = 0; i < 3; i++) this->dt[i].create (NLOCAL, NLOCAL);

    hm.hon.UOM.build_ST (1, 'T', true);

    cout << "\n===========================================================================" << endl;
    out.printrm("Force_LCAO_gamma, X Directional Derivatives of T Matrix in PW", this->dt[0]);
    cout << "\n===========================================================================" << endl;

    for (int i = 0; i < 3; i++) this->dt[i].freemem ();
    delete[] this->dt;
    */
    return;
}

// force due to the overlap matrix.
// need energy density matrix here.
void Force_LCAO_gamma::cal_foverlap(void)
{
	TITLE("Force_LCAO_gamma","cal_foverlap");
	timer::tick("Force_LCAO_gamma","cal_foverlap",'G');

	// set the force arrays.
	for(int iat=0; iat<ucell.nat; iat++)
	{
		ZEROS( foverlap[iat], 3);
	}
	for(int ipol=0; ipol<3; ++ipol)
	{
		ZEROS( this->soverlap[ipol], 3);
	}

	// set energy density matrix.
	double** edm2d = new double*[NSPIN];
	for(int is=0; is<NSPIN; ++is)
	{
		edm2d[is] = new double[ParaO.nloc];
		ZEROS(edm2d[is],ParaO.nloc);
	}

	bool with_energy = true;

	// for different spins. mohan fix bugs 2012-07-19
	for(int is=0; is<NSPIN; ++is)
	{
		this->set_EDM_gamma(edm2d[is], with_energy,is);
		//stringstream ss;
		//ss<<"EDM_"<<MY_RANK+1;
		//ofstream ofs;
		//ofs.open(ss.str().c_str());
		//ofs<<"ParaO.nrow: "<<ParaO.nrow<<"\t"<<"ParaO.ncol: "<<ParaO.ncol<<endl;
		//for(int i=0; i<ParaO.nloc; i++)
		//{
		//	if(i != 0 && i%10 == 0) ofs<<endl;
		//	ofs<<edm2d[is][i]<<"\t";
		//}
	}

	//summation \sum_{i,j} E(i,j)*dS(i,j)
	//BEGIN CALCULATION OF FORCE OF EACH ATOM


	for(int i=0; i<NLOCAL; i++)
	{
		const int iat = ucell.iwt2iat[i];
		for(int j=0; j<NLOCAL; j++)
		{
			const int mu = ParaO.trace_loc_row[j];
			const int nu = ParaO.trace_loc_col[i];

			if (mu >= 0 && nu >= 0 )
			{
				const int index = mu * ParaO.ncol + nu;

				//================================================================
				// here is the normal order, the force of each atom is calculated
				// according to each 'column' in DSloc_x, y, z
				// because the DSloc_x,y,z are anti-symmetry matrix.
				//================================================================

				double sum = 0.0;
				for(int is=0; is<NSPIN; ++is)
				{
					sum += edm2d[is][index];
				}
				sum *= 2.0;

				this->foverlap[iat][0] += sum * LM.DSloc_x[index];
				this->foverlap[iat][1] += sum * LM.DSloc_y[index];
				this->foverlap[iat][2] += sum * LM.DSloc_z[index];
				if(STRESS)
				{
					this->soverlap[0][0] += sum/2.0 * LM.DSloc_11[index];
					this->soverlap[0][1] += sum/2.0 * LM.DSloc_12[index];
					this->soverlap[0][2] += sum/2.0 * LM.DSloc_13[index];
					this->soverlap[1][1] += sum/2.0 * LM.DSloc_22[index];
					this->soverlap[1][2] += sum/2.0 * LM.DSloc_23[index];
					this->soverlap[2][2] += sum/2.0 * LM.DSloc_33[index];	
				}
			}
		}
	}

	if(STRESS)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				if(i<j) this->soverlap[j][i] = this->soverlap[i][j];
			}
		}
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				this->soverlap[i][j] *=  ucell.lat0 / ucell.omega;
			}
		}
	}
	for(int is=0; is<NSPIN; ++is)
	{
		delete[] edm2d[is];
	}
	delete[] edm2d;
	timer::tick("Force_LCAO_gamma","cal_foverlap",'G');
	return;
}

void Force_LCAO_gamma::cal_ftvnl_dphi(double** dm2d)
{	
	TITLE("Force_LCAO_gamma","cal_ftvnl_dphi");
	timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi",'G');
	if(STRESS)for(int ipol=0; ipol<3; ++ipol)
	{
		ZEROS( this->stvnl_dphi[ipol], 3);
	}
	for(int i=0; i<NLOCAL; i++)
	{
		const int iat = ucell.iwt2iat[i];
		for(int j=0; j<NLOCAL; j++)
		{
			const int mu = ParaO.trace_loc_row[j];
			const int nu = ParaO.trace_loc_col[i];

			if (mu >= 0 && nu >= 0 )
			{
				const int index = mu * ParaO.ncol + nu;
				//contribution from deriv of AO's in T+VNL term
				
				double sum = 0.0;
				for(int is=0; is<NSPIN; ++is)
				{
					sum += dm2d[is][index];
				}
				sum *= 2.0;

				this->ftvnl_dphi[iat][0] += sum * LM.DHloc_fixed_x[index];
				this->ftvnl_dphi[iat][1] += sum * LM.DHloc_fixed_y[index];
				this->ftvnl_dphi[iat][2] += sum * LM.DHloc_fixed_z[index];
				if(STRESS)
				{
					this->stvnl_dphi[0][0] += sum/2.0 * LM.DHloc_fixed_11[index];
					this->stvnl_dphi[0][1] += sum/2.0 * LM.DHloc_fixed_12[index];
					this->stvnl_dphi[0][2] += sum/2.0 * LM.DHloc_fixed_13[index];
					this->stvnl_dphi[1][1] += sum/2.0 * LM.DHloc_fixed_22[index];
					this->stvnl_dphi[1][2] += sum/2.0 * LM.DHloc_fixed_23[index];
					this->stvnl_dphi[2][2] += sum/2.0 * LM.DHloc_fixed_33[index];	
				}
			}
		}
	}
	if(STRESS){
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				if(i<j) this->stvnl_dphi[j][i] = this->stvnl_dphi[i][j];
			}
		}
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				this->stvnl_dphi[i][j] *=  ucell.lat0 / ucell.omega;
			}
		}
	}
	timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi",'G');
	return;
}

void Force_LCAO_gamma::test_gamma(double* mm, const string &name)
{
	cout << "\n PRINT " << name << endl;
	cout << setprecision(6) << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			if( abs(mm[i*NLOCAL+j])>1.0e-5)
			{
				cout << setw(12) << mm[i*NLOCAL+j];
			}
			else
			{
				cout << setw(12) << "0";
			}
		}
		cout << endl;
	}
	return;
}

void Force_LCAO_gamma::cal_fvna(LCAO_Matrix &LM)
{
	TITLE("Force_LCAO_gamma","cal_fvna");
	timer::tick("Force_LCAO_gamma","cal_fvna",'H');
	int istep=1;	
	bool delta_vh = 0;
	bool vna = 1; // tmp by mohan

	int dense = VNA;
	OUT(ofs_running,"dense grid for VNA force",dense);
	Grid_Technique gtf;
	gtf.set_pbc_grid(
	pw.ncx*dense,pw.ncy*dense,pw.ncz*dense,
	pw.bx*dense,pw.by*dense,pw.bz*dense,
	pw.nbx,pw.nby,pw.nbz,
	pw.nbxx,pw.nbzp_start,pw.nbzp,
	vna);

//--------------------------------------------------------
// The neutral potential can be only used when vna = 1,
// then to be used to check if the real space neutral
// potential is correct
//--------------------------------------------------------
//	pot.init_pot(istep, delta_vh, vna);
//	for(int ir=0; ir<pw.nrxx; ir++)
//	{
//		pot.vrs1[ir] = pot.vrs(CURRENT_SPIN, ir);
//	}

	UHM.GG.cal_force_vna(pot.vrs1, gtf, LM);
	timer::tick("Force_LCAO_gamma","cal_fvna",'H');
	return;
}

void Force_LCAO_gamma::cal_fvl_dphi(double** dm2d)
{	
	TITLE("Force_LCAO_gamma","cal_ftvnl_dphi");
	timer::tick("Force_LCAO_gamma","cal_fvl_dphi",'G');

	ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_z, ParaO.nloc);
	if(STRESS)
	{
		ZEROS (LM.DHloc_fixed_11, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_12, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_13, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_22, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_23, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_33, ParaO.nloc);
		for(int ipol=0; ipol<3; ++ipol)
		{
			ZEROS( this->svl_dphi[ipol], 3);
		}
	}

	//xiaohui add 'OUT_LEVEL', 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"VNA",VNA);

	//===================================
	// the neutral potential part.	
	//===================================
	if(VNA)
	{
		this->cal_fvna(LM);
	}

	double* tmpDHx = new double[ParaO.nloc];
	double* tmpDHy = new double[ParaO.nloc];
	double* tmpDHz = new double[ParaO.nloc];
	ZEROS( tmpDHx, ParaO.nloc );
	ZEROS( tmpDHy, ParaO.nloc );
	ZEROS( tmpDHz, ParaO.nloc );
	for(int i=0; i<ParaO.nloc; ++i)
	{
		tmpDHx[i] = LM.DHloc_fixed_x[i];
		tmpDHy[i] = LM.DHloc_fixed_y[i];
		tmpDHz[i] = LM.DHloc_fixed_z[i];
		//cout << "  LM.DHloc_fixed_x=" <<  LM.DHloc_fixed_x[i] << endl;
		//cout << "  LM.DHloc_fixed_y=" <<  LM.DHloc_fixed_y[i] << endl;
		//cout << "  LM.DHloc_fixed_z=" <<  LM.DHloc_fixed_z[i] << endl;
	}

    //calculate dVL
    //calculate <dphi | VL | phi>

	int istep = 1;
	if(VNA)
	{
		// calculate the (delta_Vh + Vxc)
		bool delta_vh = 1;
		bool vna = 0;
		// istep must set to 1, other wise the charge
		// density may be changed!
		pot.init_pot(istep, delta_vh, vna);
	}
	else
	{
		pot.init_pot(istep);
	}

	for(int is=0; is<NSPIN; ++is)
	{
		CURRENT_SPIN = is;

		ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
    	ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
    	ZEROS (LM.DHloc_fixed_z, ParaO.nloc);

		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			pot.vrs1[ir] = pot.vrs(CURRENT_SPIN, ir);
		}

		//  should not be set zero if VNA is used.
		//	ZEROS(LM.DHloc_fixed_x,ParaO.nloc);
		//	ZEROS(LM.DHloc_fixed_y,ParaO.nloc);
		//	ZEROS(LM.DHloc_fixed_z,ParaO.nloc);
		UHM.GG.cal_force(pot.vrs1);


		for(int i=0; i<NLOCAL; i++)
		{
			const int iat = ucell.iwt2iat[i];
			for(int j=0; j<NLOCAL; j++)
			{
				const int iat2 = ucell.iwt2iat[j];
				const int mu = ParaO.trace_loc_row[j];
				const int nu = ParaO.trace_loc_col[i];
				if (mu >= 0 && nu >= 0 )
				{
					const int index = mu * ParaO.ncol + nu;
					//contribution from deriv of AO's in T+VNL term

					double dm2d2 = 2.0 * dm2d[is][index];

					this->fvl_dphi[iat][0] -= dm2d2 * ( LM.DHloc_fixed_x[index] + tmpDHx[index] );
					this->fvl_dphi[iat][1] -= dm2d2 * ( LM.DHloc_fixed_y[index] + tmpDHy[index] );
					this->fvl_dphi[iat][2] -= dm2d2 * ( LM.DHloc_fixed_z[index] + tmpDHz[index] );
					if(STRESS)
					{
						this->svl_dphi[0][0] += dm2d[is][index] * LM.DHloc_fixed_11[index];
						this->svl_dphi[0][1] += dm2d[is][index] * LM.DHloc_fixed_12[index];
						this->svl_dphi[0][2] += dm2d[is][index] * LM.DHloc_fixed_13[index];
						this->svl_dphi[1][1] += dm2d[is][index] * LM.DHloc_fixed_22[index];
						this->svl_dphi[1][2] += dm2d[is][index] * LM.DHloc_fixed_23[index];
						this->svl_dphi[2][2] += dm2d[is][index] * LM.DHloc_fixed_33[index];
					}
					//	cout << setw(5) << iat << setw(5) << iat2 
					//	<< setw(5) << mu << setw(5) << nu
					//	<< setw(15) << LM.DHloc_fixed_z[index] << endl;
				}
			}
		}

//			cout << "fvl_dphi:" << endl;
//			for(int iat=0; iat<ucell.nat; ++iat)
//			{
//				cout << setw(5) << iat << setw(15) << fvl_dphi[iat][0] 
//				<< setw(15) << fvl_dphi[iat][1]
//				<< setw(15) << fvl_dphi[iat][2] << endl;
//			}


	} // end spin
	// test mohan tmp
//	test_gamma(LM.DHloc_fixed_x,"LM.DHloc_fixed_x");

	if(STRESS)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				if(i<j) this->svl_dphi[j][i] = this->svl_dphi[i][j];
			}
		}
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				this->svl_dphi[i][j] /= ucell.omega;
			}
		}
	}

	delete[] tmpDHx;
	delete[] tmpDHy;
	delete[] tmpDHz;


	timer::tick("Force_LCAO_gamma","cal_fvl_dphi",'G');
	return;
}


void Force_LCAO_gamma::cal_fvnl_dbeta(double** dm2d)
{
	TITLE("Force_LCAO_gamma","cal_fvnl_dbeta");
	timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta");
	if(STRESS)for(int ipol=0; ipol<3; ++ipol)
	{
		ZEROS( this->svnl_dbeta[ipol], 3);
	}
	for(int iat=0; iat<ucell.nat; iat++)
	{
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
		const Vector3<double> tau0 = ucell.atoms[it].tau[ia];
        //find ajacent atom of atom ia
        //GridD.Find_atom( ucell.atoms[it].tau[ia] );
        GridD.Find_atom( ucell.atoms[it].tau[ia] ,it, ia);

        //FOLLOWING ARE CONTRIBUTIONS FROM
        //VNL DUE TO PROJECTOR'S DISPLACEMENT
        for (int ad1 =0 ; ad1 < GridD.getAdjacentNum()+1; ad1++)
        {
			const int T1 = GridD.getType (ad1);
			const Atom* atom1 = &ucell.atoms[T1];
			const int I1 = GridD.getNatom (ad1);
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            const Vector3<double> tau1 = GridD.getAdjacentTau (ad1);

            for (int ad2 =0 ; ad2 < GridD.getAdjacentNum()+1; ad2++)
            {
				const int T2 = GridD.getType (ad2);
				const Atom* atom2 = &ucell.atoms[T2];
				const int I2 = GridD.getNatom (ad2);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                const Vector3<double> tau2 = GridD.getAdjacentTau (ad2);

                const double Rcut_Beta = ORB.Beta[it].get_rcut_max();
                const double Rcut_AO1 = ORB.Phi[T1].getRcut();
                const double Rcut_AO2 = ORB.Phi[T2].getRcut();

                const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                const double dist2 = (tau2-tau0).norm() * ucell.lat0;
                double r0[3],r1[3];
		if(STRESS)
                {
                    r1[0] = ( tau1.x - tau0.x) ;
                    r1[1] = ( tau1.y - tau0.y) ;
                    r1[2] = ( tau1.z - tau0.z) ;
                    r0[0] = ( tau2.x - tau0.x) ;
                    r0[1] = ( tau2.y - tau0.y) ;
                    r0[2] = ( tau2.z - tau0.z) ;
                }

                if (dist1 > Rcut_Beta + Rcut_AO1
                        || dist2 > Rcut_Beta + Rcut_AO2)
                {
                    continue;
                }

				for (int jj = 0; jj < ucell.atoms[T1].nw; jj++)
                {
                    const int iw1_all = start1 + jj;
					const int mu = ParaO.trace_loc_row[iw1_all];
					if(mu<0) continue;
                    for (int kk = 0; kk < ucell.atoms[T2].nw; kk++)
                    {
                        const int iw2_all = start2 + kk;
						const int nu = ParaO.trace_loc_col[iw2_all];
						if(nu<0) continue;
					
                        double nlm[3] = {0,0,0};
								
						UOT.snap_psibeta(
                            nlm, 1,
                            tau1, T1,
                            atom1->iw2l[jj], // L2
                            atom1->iw2m[jj], // m2
                            atom1->iw2n[jj], // N2
                            tau2, T2,
                            atom2->iw2l[kk], // L1
                            atom2->iw2m[kk], // m1
                            atom2->iw2n[kk], // n1
                            tau0, it);
						double nlm1[3] = {0,0,0};
						if(STRESS) UOT.snap_psibeta(
								nlm1, 1,
								tau2, T2,
								atom2->iw2l[kk], // L2
								atom2->iw2m[kk], // m2
								atom2->iw2n[kk], // N2
								tau1, T1,
								atom1->iw2l[jj], // L1
								atom1->iw2m[jj], // m1
								atom1->iw2n[jj], // n1
								tau0, it);

						const int index = mu * ParaO.ncol + nu;

						// dbeta is minus, that's consistent.
						// only one projector for each atom force.

						double sum = 0.0;
						for(int is=0; is<NSPIN; ++is)
						{
							sum += dm2d[is][index];
						}
						sum *= 2.0;

						this->fvnl_dbeta[iat][0] -= sum * nlm[0];
						this->fvnl_dbeta[iat][1] -= sum * nlm[1];
						this->fvnl_dbeta[iat][2] -= sum * nlm[2];

						if(STRESS) 
						{
							for(int ipol=0;ipol<3;ipol++){
								this->svnl_dbeta[0][ipol] -= sum/2.0 * (nlm[0] * r0[ipol] + nlm1[0] * r1[ipol])* -1;
								this->svnl_dbeta[1][ipol] -= sum/2.0 * (nlm[1] * r0[ipol] + nlm1[1] * r1[ipol])* -1;
								this->svnl_dbeta[2][ipol] -= sum/2.0 * (nlm[2] * r0[ipol] + nlm1[2] * r1[ipol])* -1;
							}
						}
                    }//!kk
                }//!ad2
            }//!jj
        }//!ad1
    }//!iat
	if(STRESS)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				this->svnl_dbeta[i][j] *=  ucell.lat0 / ucell.omega;
			}
		}
	}
	timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta");
	return;
}

void Force_LCAO_gamma::DerivS_PW (void)
{
    WARNING_QUIT("Force_LCAO_gamma::DerivT_PW","no use for a long time.");
    // no use for a long time.
    // mohan 2010-08-10
    /*
    //-----------------------------------------
    //TEST OUTPUT OF S Deriv Matrix in PW
    //-----------------------------------------
    //test under pw basis
    hm.hon.UOM.ForcedS = new matrix[3];
    this->ds = hm.hon.UOM.ForcedS;
    for (int i = 0; i < 3; i++) this->ds[i].create (NLOCAL, NLOCAL);

    hm.hon.UOM.build_ST (1, 'S', true);

    //    cout << "\n===========================================================================" << endl;
    //    out.printrm("Force_LCAO_gamma, X Directional Derivatives of S Matrix in PW", this->ds[2]);
    //    cout << "\n===========================================================================" << endl;

    for (int i = 0; i < 3; i++) this->ds[i].freemem ();
    delete[] this->ds;
    */
    return;
}

