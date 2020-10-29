//==========================================================
// Author:Xin Qu
// DATE : 2019-12-10
//==========================================================
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <sstream>
#include <complex>

#include "dftu.h"
#include "../input.h"
#include "../src_global/constants.h"
#include "../src_pw/global.h"
#include "global_fp.h"
#include "../src_global/global_function.h"
#include "../src_global/scalapack_connector.h"
#include "../src_global/lapack_connector.h"
#include "../src_pw/inverse_matrix.h"
#include "local_orbital_ions.h"
#include "lcao_matrix.h"
#include "../src_pw/magnetism.h"
#include "use_overlap_table.h"
#include "../src_pw/charge.h"


extern "C"
{
	void sphbsl_(int *n, double *r, double *A, double *val);
	void sphhnk_(int *n, double *r, double *A, double *val);
}


DFTU dftu;

DFTU::DFTU(){}

DFTU::~DFTU(){}


void DFTU::init()
{
    TITLE("DFTU", "init()");

	if(INPUT.dftu_type==1 && INPUT.double_counting==1) cal_type = 1;
	else if(INPUT.dftu_type==1 && INPUT.double_counting==2) cal_type = 2;
	else if(INPUT.dftu_type==2 && INPUT.double_counting==1) cal_type = 3;
	else if(INPUT.dftu_type==2 && INPUT.double_counting==2) cal_type = 4;
	else WARNING_QUIT("DFT+U", "Wrong parameter");
		
	this->EU = 0.0;

	if(FORCE)
	{
		this->force_dftu.resize(ucell.nat);
		for(int ia=0; ia<ucell.nat; ia++)
		{
			this->force_dftu.at(ia).resize(3, 0.0);
		}
	}

	if(STRESS)
	{
		this->stress_dftu.resize(3);
		for(int dim=0; dim<3; dim++)
		{
			this->stress_dftu.at(dim).resize(3, 0.0);
		}
	}
	
	if(GAMMA_ONLY_LOCAL)
	{
		this->pot_eff_gamma.resize(kv.nks);
		for(int ik=0; ik<kv.nks; ik++)
		{
			this->pot_eff_gamma.at(ik).resize(ParaO.nloc, 0.0);
		}
	}
	else
	{
		this->Sm_k.resize(kv.nks);
		this->pot_eff_k.resize(kv.nks);

		for(int ik=0; ik<kv.nks; ik++)
		{
			this->Sm_k.at(ik).resize(ParaO.nloc, ZERO);
			this->pot_eff_k.at(ik).resize(ParaO.nloc, ZERO);
		} 	
	}

	this->locale.resize(ucell.nat);
	this->locale_save.resize(ucell.nat);

	this->iatlnmipol2iwt.resize(ucell.nat);
	this->iat2it.resize(ucell.nat);
	this->iwt2it.resize(NLOCAL);
	this->iwt2iat.resize(NLOCAL);
	this->iwt2l.resize(NLOCAL);
	this->iwt2n.resize(NLOCAL);
	this->iwt2m.resize(NLOCAL);
	this->iwt2ipol.resize(NLOCAL);	
	for(int i=0; i<NLOCAL; i++)
	{
		this->iwt2it.at(i) = -1;
		this->iwt2iat.at(i) = -1;
		this->iwt2l.at(i) = -1;
		this->iwt2n.at(i) = -1;
		this->iwt2m.at(i) = -1;
		this->iwt2ipol.at(i) = -1;
	}

	for(int it=0; it<INPUT.ntype; ++it)
	{				
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			const int iat = ucell.itia2iat(it, ia);
			this->iat2it.at(iat) = it;

			locale.at(iat).resize(ucell.atoms[it].nwl+1);
			locale_save.at(iat).resize(ucell.atoms[it].nwl+1);

			for(int l=0; l<=ucell.atoms[it].nwl; l++)
			{			
				const int N = ucell.atoms[it].l_nchi[l];

				locale.at(iat).at(l).resize(N);
				locale_save.at(iat).at(l).resize(N);

				for(int n=0; n<N; n++)
				{
					locale.at(iat).at(l).at(n).resize(2);
					locale_save.at(iat).at(l).at(n).resize(2);

					locale.at(iat).at(l).at(n).at(0).create(2*l+1, 2*l+1);
					locale.at(iat).at(l).at(n).at(1).create(2*l+1, 2*l+1);

					locale_save.at(iat).at(l).at(n).at(0).create(2*l+1, 2*l+1);
					locale_save.at(iat).at(l).at(n).at(1).create(2*l+1, 2*l+1);												
				}
			}

			//initialize the arrry iatlnm2iwt[iat][l][n][m]
			this->iatlnmipol2iwt.at(iat).resize(ucell.atoms[it].nwl+1);
			for(int L=0; L<=ucell.atoms[it].nwl; L++)
			{
				this->iatlnmipol2iwt.at(iat).at(L).resize(ucell.atoms[it].l_nchi[L]);

				for(int n=0; n<ucell.atoms[it].l_nchi[L]; n++)
				{
					this->iatlnmipol2iwt.at(iat).at(L).at(n).resize(2*L+1);
					
					for(int m=0; m<2*L+1; m++)
					{
						this->iatlnmipol2iwt.at(iat).at(L).at(n).at(m).resize(NPOL);
					}
				}
			}

			for(int iw=0; iw<ucell.atoms[it].nw*NPOL; iw++)
			{
				int iw0 = iw/NPOL;
				int ipol = iw%NPOL;
				int iwt = ucell.itiaiw2iwt(it, ia, iw);
				int l = ucell.atoms[it].iw2l[iw0];
				int n = ucell.atoms[it].iw2n[iw0];
				int m = ucell.atoms[it].iw2m[iw0];
								
				this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m).at(ipol) = iwt;
				this->iwt2it.at(iwt) = it;
				this->iwt2iat.at(iwt) = iat;
				this->iwt2l.at(iwt) = l;
				this->iwt2n.at(iwt) = n;
				this->iwt2m.at(iwt) = m;
				this->iwt2ipol.at(iwt) = ipol;
			}
		}	
	}

	this->Yukawa = INPUT.yukawa_potential;
 	if(Yukawa)
 	{
		this->Fk.resize(INPUT.ntype);
		
		this->U_Yukawa.resize(INPUT.ntype);
		this->J_Yukawa.resize(INPUT.ntype);	

		for(int it=0; it<ucell.ntype; it++)
		{			
			const int NL = ucell.atoms[it].nwl + 1;

			this->Fk.at(it).resize(NL);		
			this->U_Yukawa.at(it).resize(NL);
			this->J_Yukawa.at(it).resize(NL);	

			for(int l=0; l<NL; l++)
			{
				int N = ucell.atoms[it].l_nchi[l];

				this->Fk.at(it).at(l).resize(N);
				for(int n=0; n<N; n++)
				{
					this->Fk.at(it).at(l).at(n).resize(l+1, 0.0);
				}	

				this->U_Yukawa.at(it).at(l).resize(N, 0.0);
				this->J_Yukawa.at(it).at(l).resize(N, 0.0);

				// if(l>=INPUT.orbital_corr[it] && INPUT.orbital_corr[it]!=-1)
				// {
					// this->cal_slater_Fk(l, it);
					// this->cal_slater_UJ(it, l);
				// }				
			}			 	
		}
 	}
	else
	{
		this->U = INPUT.hubbard_u;            //Hubbard Coulomb interaction parameter U(ev)
		this->J = INPUT.hund_j;               //Hund exchange parameter J(ev)
	}

	if(CALCULATION=="nscf")
	{
		stringstream sst; 
		sst << global_out_dir << "onsite.dm"; 
		this->read_occup_m( sst.str() );
		this->local_occup_bcast();		
	}
	else
	{
		if(INPUT.omc) 
		{
			stringstream sst; 
			sst << "initial_onsite.dm"; 
			this->read_occup_m( sst.str() );
			this->local_occup_bcast();
		}
	}

	//this->out_numorb();

  //ofs_running << "dftu.cpp "<< __LINE__ << endl;
    return;
}


void DFTU::cal_occup_m_k(const int iter)
{
	TITLE("DFTU", "cal_occup_m_k");

	for(int T=0; T<ucell.ntype; T++)
	{
		if(INPUT.orbital_corr[T]==-1) continue;

		//const int LC = INPUT.orbital_corr[T];
		
		for(int I=0; I<ucell.atoms[T].na; I++)
		{
		 	const int iat = ucell.itia2iat(T, I);
			
			for(int l=0; l<ucell.atoms[T].nwl+1; l++)
			{
				const int N = ucell.atoms[T].l_nchi[l];

				for(int n=0; n<N; n++)
				{						
					locale_save.at(iat).at(l).at(n).at(0) = locale.at(iat).at(l).at(n).at(0);
					locale_save.at(iat).at(l).at(n).at(1) = locale.at(iat).at(l).at(n).at(1);
				}
			}			
		}
	}

	//=================Part 1======================
	//call PBLAS routine to calculate the product of the S and density matrix
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;

	vector<vector<complex<double>>> srho(kv.nks);
	for(int ik=0; ik<kv.nks; ik++)
	{
		srho.at(ik).resize(ParaO.nloc, complex<double>(0.0, 0.0));
	}	
	
	for(int ik=0; ik<kv.nks; ik++)
	{		
		// srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
		pzgemm_(&transN, &transT,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				VECTOR_TO_PTR(Sm_k.at(ik)), &one_int, &one_int, ParaO.desc, 
				LOC.wfc_dm_2d.dm_k.at(ik).c, &one_int, &one_int, ParaO.desc,
				&beta, 
				VECTOR_TO_PTR(srho.at(ik)), &one_int, &one_int, ParaO.desc);
				
	}


	//test the sum rule of srho
	/*
	complex<double> elec_tot(0.0, 0.0);
	complex<double> Nele(0.0, 0.0);
	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				int row = ParaO.MatrixInfo.row_set[ir];
				int col = ParaO.MatrixInfo.col_set[ic];

				if(row==col) elec_tot += srho.at(ik).at(ic*ParaO.nrow + ir);
			}
		}
	}
	MPI_Allreduce(&elec_tot, &Nele, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
 	this->Nval = Nele.real();
	if(MY_RANK==0)
	{
		ofstream ofs_elec("nelec.dat", ios_base::app);
		ofs_elec << "Total number of electrons of the system " << Nele.real() << " i" << Nele.imag() << endl;
	}
	*/

	//=================Part 2======================
	//get the local occupation number matrix from the product of S and density matrix, i.e. srho。
	for(int it=0; it<ucell.ntype; it++)
	{
		const int NL = ucell.atoms[it].nwl + 1;
		const int LC = INPUT.orbital_corr[it];
		
		if(LC == -1) continue;

		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{		
			const int iat = ucell.itia2iat(it, ia);

			for(int l=0; l<NL; l++)
			{				
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[it]) continue;
				// }
				// else
				// {
					// if(l!=INPUT.orbital_corr[it]) continue;
				// }
				if(l!=INPUT.orbital_corr[it]) continue;

				const int N = ucell.atoms[it].l_nchi[l];
						
				for(int n=0; n<N; n++)
				{
				 	// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;
					// set the local occupation mumber matrix of spin up and down zeros					
					locale.at(iat).at(l).at(n).at(0).zero_out();
					locale.at(iat).at(l).at(n).at(1).zero_out();

					vector<ComplexMatrix> loc_occup_m;
					vector<ComplexMatrix> loc_occup_m_tmp;
					if(NSPIN==1)
					{
						loc_occup_m.resize(1);
						loc_occup_m_tmp.resize(1);

						loc_occup_m.at(0).create((2*l+1)*NPOL, (2*l+1)*NPOL);
						loc_occup_m_tmp.at(0).create((2*l+1)*NPOL, (2*l+1)*NPOL);
					}
					else if(NSPIN==2 || NSPIN==4)
					{
						loc_occup_m.resize(2);
						loc_occup_m_tmp.resize(2);

						loc_occup_m.at(0).create(2*l+1, 2*l+1);
						loc_occup_m_tmp.at(0).create(2*l+1, 2*l+1);

						loc_occup_m.at(1).create(2*l+1, 2*l+1);
						loc_occup_m_tmp.at(1).create(2*l+1, 2*l+1);
					}		

					//Calculate the local occupation number matrix			
					for(int m0=0; m0<2*l+1; m0++)
					{
						for(int ipol0=0; ipol0<NPOL; ipol0++)
						{
							const int iwt0 = this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m0).at(ipol0);
							const int mu = ParaO.trace_loc_row[iwt0];
							const int mu_prime = ParaO.trace_loc_col[iwt0];

							for(int m1=0; m1<2*l+1; m1++)
							{
								for(int ipol1=0; ipol1<NPOL; ipol1++)
								{									
									const int iwt1 = this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m1).at(ipol1);
									const int nu = ParaO.trace_loc_col[iwt1];
									const int nu_prime = ParaO.trace_loc_row[iwt1];

									const int irc = nu*ParaO.nrow + mu;
									const int irc_prime = mu_prime*ParaO.nrow + nu_prime;

									// const int m0_all = m0 + ipol0*(2*l+1);
									// const int m1_all = m1 + ipol1*(2*l+1);

									if( (nu>=0) && (mu>=0) )
									{	
										if(NSPIN==1 || NSPIN==2)
										{
											for(int ik=0; ik<kv.nks; ik++)
											{
												const int spin = kv.isk[ik];

												loc_occup_m_tmp.at(spin)(m0, m1) += srho.at(ik).at(irc)/4.0;
											}
										}	
										else if(NSPIN==4)//SOC
										{
											if(ipol0==ipol1)
											{
												const int spin = ipol0;
												for(int ik=0; ik<kv.nks; ik++)
												{
													loc_occup_m_tmp.at(spin)(m0, m1) += srho.at(ik).at(irc)/4.0;
												} 
											}
										}												
									}

									if( (nu_prime>=0) && (mu_prime>=0) )
									{
										if(NSPIN==1 || NSPIN==2)
										{
											for(int ik=0; ik<kv.nks; ik++)
											{
												const int spin = kv.isk[ik];

												loc_occup_m_tmp.at(spin)(m0, m1) += std::conj(srho.at(ik).at(irc_prime))/4.0;
											}
										}
										else if(NSPIN==4) //SOC
										{
											if(ipol0==ipol1)
											{
												const int spin = ipol0;
												for(int ik=0; ik<kv.nks; ik++)
												{
													loc_occup_m_tmp.at(spin)(m0, m1) += std::conj(srho.at(ik).at(irc_prime))/4.0;
												}
											}
										}
									}								
								}//ipol1										
							}//m1
						}//ipol0
					}//m0
					
					for(int m0=0; m0<2*l+1; m0++)
					{
						for(int m1=0; m1<2*l+1; m1++)
						{
							if(NSPIN==1 )
							{
								MPI_Allreduce( &loc_occup_m_tmp.at(0)(m0,m1), &loc_occup_m.at(0)(m0,m1), 1,
												MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
							}
							else if(NSPIN==2 || NSPIN==4)
							{
								MPI_Allreduce( &loc_occup_m_tmp.at(0)(m0,m1), &loc_occup_m.at(0)(m0,m1), 1,
												MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
								
								MPI_Allreduce( &loc_occup_m_tmp.at(1)(m0,m1), &loc_occup_m.at(1)(m0,m1), 1,
												MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
							}
						}
					}
				
					// for the case spin independent calculation
					switch(NSPIN)
					{
					case 1:
						for(int is=0; is<2; is++)
						{
							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int m1=0; m1<2*l+1; m1++)
								{
									locale.at(iat).at(l).at(n).at(is)(m0,m1) = 0.5*(loc_occup_m.at(0)(m0,m1).real() + loc_occup_m.at(0)(m1,m0).real());
								}			
							}	
						}
						break;

					case 2:
						for(int is=0; is<NSPIN; is++)
						{
							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int m1=0; m1<2*l+1; m1++)
								{
									locale.at(iat).at(l).at(n).at(is)(m0,m1) = loc_occup_m.at(is)(m0,m1).real() + loc_occup_m.at(is)(m1,m0).real();
								}		
							}
						}
						break;

					case 4: //SOC
						for(int m0=0; m0<2*l+1; m0++)
						{							
							for(int m1=0; m1<2*l+1; m1++)
							{								
								locale.at(iat).at(l).at(n).at(0)(m0, m1) = loc_occup_m.at(0)(m0, m1).real() + loc_occup_m.at(0)(m1, m0).real();	
								locale.at(iat).at(l).at(n).at(1)(m0, m1) = loc_occup_m.at(1)(m0, m1).real() + loc_occup_m.at(1)(m1, m0).real();																		
							}										
						}
						break;

					default:
						cout << "Not supported NSPIN parameter" << endl;
						exit(0);			
					}

						//test
						/*
						if(MY_RANK==0)
						{
							if(NSPIN==4)
							{
								ofstream of_soc("soc_locmat.dat", ios_base::app);

								for(int m0=0; m0<2*l+1; m0++)
								{
									for(int ipol0=0; ipol0<NPOL; ipol0++)
									{
										const int m0_all = m0 + (2*l+1)*ipol0;

										for(int m1=0; m1<2*l+1; m1++)
										{
											for(int ipol1=0; ipol1<NPOL; ipol1++)
											{
												const int m1_all = m1 + (2*l+1)*ipol1;

												complex<double> nmm = 0.5*(loc_occup_m.at(0)(m0_all, m1_all) + loc_occup_m.at(0)(m1_all, m0_all));

												of_soc << "(" << fixed << setw(8) << setprecision(4) << nmm.real() << " i" 
												<< fixed << setw(8) << setprecision(4) << nmm.imag() << ")    ";
											}							
										}
										of_soc << endl;
									}	
								}	
								of_soc << "TWO_FERMI  " << TWO_EFERMI << endl;
								of_soc << endl;
								of_soc << endl;
							}
							else if(NSPIN==2)
							{
								ofstream of_soc("nonsoc_locmat.dat", ios_base::app);

								for(int m0=0; m0<2*l+1; m0++)
								{
									for(int ipol0=0; ipol0<NPOL; ipol0++)
									{
										const int m0_all = m0 + (2*l+1)*ipol0;

										for(int m1=0; m1<2*l+1; m1++)
										{
											for(int ipol1=0; ipol1<NPOL; ipol1++)
											{
												const int m1_all = m1 + (2*l+1)*ipol1;

												complex<double> nmm = 0.5*(loc_occup_m.at(0)(m0_all, m1_all) + loc_occup_m.at(0)(m1_all, m0_all));

												of_soc << "(" << fixed << setw(8) << setprecision(4) << nmm.real() << " i" 
												<< fixed << setw(8) << setprecision(4) << nmm.imag() << ")    ";
											}							
										}
										of_soc << endl;
									}	
								}
								of_soc << "TWO_FERMI  " << TWO_EFERMI << endl;	
								of_soc << endl;
								of_soc << endl;
							}
						}
						*/

				}//end n

				// this->print(it, iat, l, N, iter);
			}//end l
		}//end ia
	}//end it

	//ofs_running << "dftu.cpp "<< __LINE__  << endl;
	return;
}


void DFTU::cal_occup_m_gamma(const int iter)
{
	TITLE("DFTU", "cal_occup_m_gamma");	

	for(int T=0; T<ucell.ntype; T++)
	{
		if(INPUT.orbital_corr[T]==-1) continue;
		const int L = INPUT.orbital_corr[T];
		
		for(int I=0; I<ucell.atoms[T].na; I++)
		{
			const int iat = ucell.itia2iat(T, I);
			
			for(int l=L; l<ucell.atoms[T].nwl+1; l++)
			{
				const int N = ucell.atoms[T].l_nchi[l];

				for(int n=0; n<N; n++)
				{									
					locale_save.at(iat).at(l).at(n).at(0) = locale.at(iat).at(l).at(n).at(0);
					locale_save.at(iat).at(l).at(n).at(1) = locale.at(iat).at(l).at(n).at(1);					
				}
			}			
		}
	}
	
	//=================Part 1======================
	//call PBLAS routine to calculate the product of the S and density matrix
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;

	vector<vector<double>> srho(NSPIN);
	for(int is=0; is<NSPIN; is++)
	{
		srho.at(is).resize(ParaO.nloc);
		ZEROS(VECTOR_TO_PTR(srho.at(is)), ParaO.nloc);
	
		// srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_gamma(iw,nu)
		pdgemm_(&transN, &transT,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				LM.Sloc, &one_int, &one_int, ParaO.desc, 
				LOC.wfc_dm_2d.dm_gamma.at(is).c, &one_int, &one_int, ParaO.desc,
				&beta,
				VECTOR_TO_PTR(srho.at(is)), &one_int, &one_int, ParaO.desc);
	}

	//test the sum rule of srho
	/*
	double elec_tot = 0.0;
	double Nele = 0.0; 
	for(int is=0; is<NSPIN; is++)
	{
		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				int row = ParaO.MatrixInfo.row_set[ir];
				int col = ParaO.MatrixInfo.col_set[ic];

				if(row==col) elec_tot += srho.at(is).at(ic*ParaO.nrow + ir);
			}
		}
	}	
	MPI_Allreduce(&elec_tot, &Nele, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if(MY_RANK==0)
	{
		ofstream ofs_elec("nelec.dat", ios_base::app);
		ofs_elec << "Total number of electrons of the system " << Nele << endl;
	}
	*/

	//=================Part 2======================
	//get the local occupation number matrix from the product of S and density matrix, i.e. srho。
	for(int it=0; it<ucell.ntype; it++)
	{
		const int NL = ucell.atoms[it].nwl + 1;
		if(INPUT.orbital_corr[it] == -1) continue;
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			const int iat = ucell.itia2iat(it, ia);

			for(int l=0; l<NL; l++)
			{				
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[it]) continue;
				// }
				// else
				// {
					// if(l!=INPUT.orbital_corr[it]) continue;
				// }
				if(l!=INPUT.orbital_corr[it]) continue;

				const int N = ucell.atoms[it].l_nchi[l];

				for(int n=0; n<N; n++)
				{
				 	// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;
					// set the local occupation mumber matrix of spin up and down zeros					
					locale.at(iat).at(l).at(n).at(0).zero_out();
					locale.at(iat).at(l).at(n).at(1).zero_out();					

					vector<matrix> loc_occup_m(NSPIN);
					vector<matrix> loc_occup_m_tmp(NSPIN);

					for(int is=0; is<NSPIN; is++)
					{
						loc_occup_m.at(is).create(2*l+1, 2*l+1);
						loc_occup_m_tmp.at(is).create(2*l+1, 2*l+1);
					}

					//Calculate the local occupation number matrix			
					for(int m0=0; m0<2*l+1; m0++)
					{	
						for(int ipol0=0; ipol0<NPOL; ipol0++)
						{
							const int iwt0 = this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m0).at(ipol0);
							const int mu = ParaO.trace_loc_row[iwt0];
							const int mu_prime = ParaO.trace_loc_col[iwt0];

							for(int m1=0; m1<2*l+1; m1++)
							{	
								for(int ipol1=0; ipol1<NPOL; ipol1++)
								{											
									const int iwt1 = this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m1).at(ipol1);
									const int nu = ParaO.trace_loc_col[iwt1];
									const int nu_prime = ParaO.trace_loc_row[iwt1];

									const int irc = nu*ParaO.nrow + mu;
									const int irc_prime = mu_prime*ParaO.nrow + nu_prime;

									for(int ik=0; ik<kv.nks; ik++)
									{
										int spin = kv.isk[ik];

										if( (nu>=0) && (mu>=0) )
										{																																																
											int m0_all = m0 + (2*l+1)*ipol0;
											int m1_all = m0 + (2*l+1)*ipol1;

											loc_occup_m_tmp.at(spin)(m0,m1) += srho.at(spin).at(irc)/4.0;														
										}

										if( (nu_prime>=0) && (mu_prime>=0) )
										{
											int m0_all = m0 + (2*l+1)*ipol0;
											int m1_all = m0 + (2*l+1)*ipol1;
											
											loc_occup_m_tmp.at(spin)(m0,m1) += srho.at(spin).at(irc_prime)/4.0;
										}
									}	
								}			
							}
						}
					}

					for(int m0=0; m0<2*l+1; m0++)
					{
						for(int m1=0; m1<2*l+1; m1++)
						{
							for(int is=0; is<NSPIN; is++)
							{
								MPI_Allreduce( &loc_occup_m_tmp.at(is)(m0,m1), &loc_occup_m.at(0)(m0,m1), 1,
												MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
							}
							
						}
					}

					// for the case spin independent calculation
					switch(NSPIN)
					{
					case 1:	
						for(int is=0; is<2; is++)
						{
							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int m1=0; m1<2*l+1; m1++)
								{
									locale.at(iat).at(l).at(n).at(is)(m0,m1) = 0.5*(loc_occup_m.at(0)(m0,m1) + loc_occup_m.at(0)(m1,m0));
								}		
							}	
						}	
						break;

					case 2:
						for(int is=0; is<NSPIN; is++)
						{
							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int m1=0; m1<2*l+1; m1++)
								{
									locale.at(iat).at(l).at(n).at(is)(m0,m1) = loc_occup_m.at(0)(m0,m1) + loc_occup_m.at(0)(m1,m0);
								}					
							}
						}
						break;

					default:
						cout << "Not supported NSPIN parameter" << endl;
						exit(0);			
					}

				}//end for(n)

				//this->print(it, iat, 2*L+1, N);
			}
		}
	}
	
	//ofs_running << "dftu.cpp "<< __LINE__  << endl;
	return;
}


void DFTU::write_occup_m(const string &fn)
{
	TITLE("DFTU", "write_occup_m");

	if(MY_RANK!=0) return;

	ofstream ofdftu;
	ofdftu.open(fn.c_str());
    if (!ofdftu)
    {
    	cout << "DFTU::write_occup_m. Can't create file onsite.dm!" << endl;
		exit(0);
    }

	for(int T=0; T<ucell.ntype; T++)
	{
		if(INPUT.orbital_corr[T]==-1) continue;
		const int NL = ucell.atoms[T].nwl+1;
		const int LC = INPUT.orbital_corr[T];

		for(int I=0; I<ucell.atoms[T].na; I++)
		{
			const int iat = ucell.itia2iat(T, I);
			ofdftu << "atoms" << "  " << iat << endl;

			for(int l=0; l<NL; l++)
			{
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=INPUT.orbital_corr[T]) continue;
				// }
				if(l!=INPUT.orbital_corr[T]) continue;
				
				const int N = ucell.atoms[T].l_nchi[l];
				ofdftu << "L" << "  " << l << endl;

				for(int n=0; n<N; n++)		
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

					ofdftu << "zeta" << "  " << n << endl;
					
					for(int is=0; is<2; is++)
					{
						ofdftu << "spin" << "  " << is << endl;
						for(int m0=0; m0<2*l+1; m0++)
						{
							for(int m1=0; m1<2*l+1; m1++)
							{
								ofdftu << fixed << setw(12) << setprecision(8) << locale.at(iat).at(l).at(n).at(is)(m0, m1);
							}
							ofdftu << endl;
						}
					}					
				}//n
			}//l
		}//I
	}//T

	return;
}


void DFTU::read_occup_m(const string &fn)
{
	TITLE("DFTU", "read_occup_m");

	if(MY_RANK!=0) return;

	ifstream ifdftu(fn.c_str(), ios::in);	

	if (!ifdftu) 
	{
		if(CALCULATION=="nscf") 
		{
			cout << "DFTU::read_occup_m. Can not find the file oneite.dm . Please do scf calculation first" << endl;
		}
		else
		{
			if(INPUT.omc) 
			{
				cout << "DFTU::read_occup_m. Can not find the file intial_oneite.dm . Please check your intial_onsite.dm" << endl;
			}
		}
		exit(0);
	}

	ifdftu.clear();
    ifdftu.seekg(0);

	char word[10];

	int T, iat, spin, L, zeta;

	ifdftu.rdstate();

    while(ifdftu.good())
    {
        ifdftu >> word;
		if(ifdftu.eof()) break;
		
		if(strcmp("atoms", word) == 0)
        {
            ifdftu >> iat;
			ifdftu.ignore(150, '\n');

			T= this->iat2it.at(iat);
			const int NL = ucell.atoms[T].nwl + 1;
			const int LC = INPUT.orbital_corr[T];
	
			for(int l=0; l<NL; l++)			
			{
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=INPUT.orbital_corr[T]) continue;
				// }
				if(l!=INPUT.orbital_corr[T]) continue;

				ifdftu >> word;

				if (strcmp("L", word) == 0)
				{
					ifdftu >> L;
					ifdftu.ignore(150, '\n');

					const int N = ucell.atoms[T].l_nchi[L];
					for(int n=0; n<N; n++)
					{
						// if(!Yukawa && n!=0) continue;
						if(n!=0) continue;

						ifdftu >> word;
						if(strcmp("zeta", word) == 0)						
						{														
							ifdftu >> zeta;
							ifdftu.ignore(150, '\n');
							
							for(int is=0; is<2; is++)
							{
								ifdftu >> word;
								if(strcmp("spin", word) == 0)
								{								
									ifdftu >> spin;
									ifdftu.ignore(150, '\n');

									double value = 0.0;
									for(int m0=0; m0<2*L+1; m0++)
									{
										for(int m1=0; m1<2*L+1; m1++)
										{							
											ifdftu >> value;
											locale.at(iat).at(L).at(zeta).at(spin)(m0, m1) = value;						 	
										}
										ifdftu.ignore(150, '\n');											
									}
								}
								else
								{								
									cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << endl;
									exit(0);
								}	
							}																																		
						}
						else
						{						
							cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << endl;
							exit(0);
						}
					}
				}
				else
				{				
					cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << endl;
					exit(0);
				}								
			}			
        }
		else
		{		
			cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << endl;
			exit(0);
		}
		
		ifdftu.rdstate();

		if (ifdftu.eof() != 0)
  		{
			break;
  		}
	}
   
	return;
}


void DFTU::local_occup_bcast()
{
	TITLE("DFTU", "local_occup_bcast");

	for(int T=0; T<ucell.ntype; T++)
	{
		if(INPUT.orbital_corr[T]==-1) continue;

		for(int I=0; I<ucell.atoms[T].na; I++)
		{
			const int iat = ucell.itia2iat(T,I);
			const int L = INPUT.orbital_corr[T];

			for(int l=0; l<=ucell.atoms[T].nwl; l++)
			{
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=INPUT.orbital_corr[T]) continue;
				// }
				if(l!=INPUT.orbital_corr[T]) continue;

				for(int n=0; n<ucell.atoms[T].l_nchi[l]; n++)
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;
										
					for(int spin=0; spin<2; spin++)
					{
					   for(int m0=0; m0<2*l+1; m0++)
					   {
						   for(int m1=0; m1<2*l+1; m1++)
						   {
								MPI_Bcast(&locale.at(iat).at(l).at(n).at(spin)(m0, m1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
						   }
					   }
					}					
				}
			}
		}
	}	    
	return;
}


void DFTU::cal_energy_correction(const int istep)
{
	TITLE("DFTU", "cal_energy_correction");

 	if((CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax") && (!INPUT.omc) && istep==0 && this->iter_dftu==1) return;

	this->EU = 0.0;
	double EU_dc = 0.0;
	
	for(int T=0; T<ucell.ntype; T++)
	{
		const int NL = ucell.atoms[T].nwl + 1;
		const int LC = INPUT.orbital_corr[T]; 
		for(int I=0; I<ucell.atoms[T].na; I++)
		{
			if(LC == -1) continue;

			const int iat = ucell.itia2iat(T, I);
			const int L = INPUT.orbital_corr[T];

			for(int l=0; l<NL; l++)
			{
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[T]) continue;
				// }	
				// else if(l!=INPUT.orbital_corr[T]) continue;
				if(l!=INPUT.orbital_corr[T]) continue;

				const int N = ucell.atoms[T].l_nchi[l];

				const int m_tot = 2*l+1;

				//part 1: calculate the DFT+U energy correction
				for(int n=0; n<N; n++)
				{
				 	// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

					if(cal_type==1) //rotationally invarient formalism and FLL double counting; Not available at present
					{

					}
					else if(cal_type==2) //rotationally invarient formalism and AMF double counting; Not available at present
					{

					}
					else if(cal_type==3) //simplified formalism and FLL double counting
					{	
						for(int spin=0; spin<2; spin++)
						{
							double nm_trace = 0.0;
							double nm2_trace = 0.0;

							for(int m0=0; m0<2*l+1; m0++)
							{
							 	nm_trace += this->locale.at(iat).at(l).at(n).at(spin)(m0, m0);
								for(int m1=0; m1<2*l+1; m1++)
								{
									nm2_trace += this->locale.at(iat).at(l).at(n).at(spin)(m0,m1)*this->locale.at(iat).at(l).at(n).at(spin)(m1,m0);
								}
							}
							if(Yukawa) this->EU += 0.5*(this->U_Yukawa.at(T).at(l).at(n) - this->J_Yukawa.at(T).at(l).at(n))*(nm_trace - nm2_trace);
							else this->EU += 0.5*(this->U[T] - this->J[T])*(nm_trace - nm2_trace);	
						}
						
					}
					else if(cal_type==4) //simplified formalism and AMF double counting; ot available at present
					{

					}

					//calculate the double counting term included in eband													
					for(int m1=0; m1<2*l+1; m1++)
					{						
						for(int m2=0; m2<2*l+1; m2++)
						{								
							for(int is=0; is<2; is++)
							{
								double VU = 0.0;
								VU = get_onebody_eff_pot(T, iat, l, n, is, m1, m2, cal_type, 0);
								EU_dc += VU*this->locale.at(iat).at(l).at(n).at(is)(m1,m2);
							}					
						}						
					}						
				}//end n
			}//end L
		}//end I
	}//end T
	
	/*
	if(MY_RANK==0)
	{
		ofstream of_eu("energy_correction", ios_base::app);
		double e=EU*Ry_to_eV;
		of_eu << "ITERATION STEP " << setw(3) << this->iter_dftu <<"      " << setw(15) << e << "ev" << endl;
	}
	*/
	
	//substract the double counting EU_dc included in band energy eband
	this->EU -= EU_dc;

	/*
	if(MY_RANK==0)
	{
		ofstream of_en("energy_correction.dat",ios_base::app);
		double val = this->EU*Ry_to_eV;

		of_en << fixed << setprecision(8) << val << "eV" << endl;
	}
	*/

	return;
}


void DFTU::cal_eff_pot_mat(const int ik, const int istep)
{
	TITLE("DFTU", "cal_eff_pot_mat");

 	if((CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax") && (!INPUT.omc) && istep==0 && this->iter_dftu==1) return;
	
	if(GAMMA_ONLY_LOCAL) ZEROS(VECTOR_TO_PTR(this->pot_eff_gamma.at(kv.isk[ik])), ParaO.nloc);
	else ZEROS(VECTOR_TO_PTR(this->pot_eff_k.at(ik)), ParaO.nloc);

	//ofs_running << "dftu.cpp "<< __LINE__  << endl;
	//=============================================================
	//   PART2: call pblas to calculate effective potential matrix
	//=============================================================
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;

	if(GAMMA_ONLY_LOCAL)
	{
		int spin = kv.isk[ik];
		vector<double> VU(ParaO.nloc, 0.0);

		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int T1 = this->iwt2it.at(iwt1);
			const int iat1 = this->iwt2iat.at(iwt1);
			const int L1 = this->iwt2l.at(iwt1);
			const int n1 = this->iwt2n.at(iwt1);
			const int m1 = this->iwt2m.at(iwt1);
			const int ipol1 = this->iwt2ipol.at(iwt1);

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int T2 = this->iwt2it.at(iwt2);
				const int iat2 = this->iwt2iat.at(iwt2);
				const int L2 = this->iwt2l.at(iwt2);
				const int n2 = this->iwt2n.at(iwt2);
				const int m2 = this->iwt2m.at(iwt2);
				const int ipol2 = this->iwt2ipol.at(iwt2);

				int irc = ic*ParaO.nrow + ir;			

				if(INPUT.orbital_corr[T1]==-1 || INPUT.orbital_corr[T2]==-1) continue;
				if(iat1!=iat2) continue;			
				// if(Yukawa)
				// {
					// if(L1<INPUT.orbital_corr[T1] || L2<INPUT.orbital_corr[T2]) continue;
				// }
				// else
				// {
					if(L1!=INPUT.orbital_corr[T1] || L2!=INPUT.orbital_corr[T2] || n1!=0 || n2!=0) continue;
				// }
				if(L1!=L2 || n1!=n2) continue;

				double val = get_onebody_eff_pot(T1, iat1, L1, n1, spin, m1, m2, cal_type, 1);

				VU.at(irc) = val;	
			}
		}
		
		vector<double> potm_tmp(ParaO.nloc, 0.0);

		// The first term
		pdgemm_(&transN, &transN,
			&NLOCAL, &NLOCAL, &NLOCAL,
			&alpha, 
			VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc, 
			LM.Sloc, &one_int, &one_int, ParaO.desc,
			&beta,
			VECTOR_TO_PTR(potm_tmp), &one_int, &one_int, ParaO.desc);
		
		for(int irc=0; irc<ParaO.nloc; irc++)
		{
			this->pot_eff_gamma.at(spin).at(irc) += 0.5*potm_tmp.at(irc);
		}

		// The second term
		ZEROS(VECTOR_TO_PTR(potm_tmp), ParaO.nloc);

		pdgemm_(&transN, &transN,
			&NLOCAL, &NLOCAL, &NLOCAL,
			&alpha, 
			LM.Sloc, &one_int, &one_int, ParaO.desc, 
			VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc,
			&beta,
			VECTOR_TO_PTR(potm_tmp), &one_int, &one_int, ParaO.desc);
		
		for(int irc=0; irc<ParaO.nloc; irc++)
		{
			this->pot_eff_gamma.at(spin).at(irc) += 0.5*potm_tmp.at(irc);
		}
	}
	else
	{
		vector<complex<double>> VU(ParaO.nloc, complex<double>(0.0, 0.0));

		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int T1 = this->iwt2it.at(iwt1);
			const int iat1 = this->iwt2iat.at(iwt1);
			const int L1 = this->iwt2l.at(iwt1);
			const int n1 = this->iwt2n.at(iwt1);
			const int m1 = this->iwt2m.at(iwt1);
			const int ipol1 = this->iwt2ipol.at(iwt1);

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int T2 = this->iwt2it.at(iwt2);
				const int iat2 = this->iwt2iat.at(iwt2);
				const int L2 = this->iwt2l.at(iwt2);
				const int n2 = this->iwt2n.at(iwt2);
				const int m2 = this->iwt2m.at(iwt2);
				const int ipol2 = this->iwt2ipol.at(iwt2);

				int irc = ic*ParaO.nrow + ir;			

				if(INPUT.orbital_corr[T1]==-1 || INPUT.orbital_corr[T2]==-1) continue;
				if(iat1!=iat2) continue;			
				// if(Yukawa)
				// {
					// if(L1<INPUT.orbital_corr[T1] || L2<INPUT.orbital_corr[T2]) continue;
				// }
				// else
				// {
					if(L1!=INPUT.orbital_corr[T1] || L2!=INPUT.orbital_corr[T2] || n1!=0 || n2!=0) continue;
				// }
				if(L1!=L2 || n1!=n2) continue;

				// if(m1==m2 && iwt1==iwt2) delta.at(irc) = 1.0;

				if(NSPIN==1 || NSPIN==2)
				{
					int spin = kv.isk[ik];
					double val = get_onebody_eff_pot(T1, iat1, L1, n1, spin, m1, m2, cal_type, 1);
					VU.at(irc) = complex<double>(val, 0.0);	
				}
				else if(NSPIN==4) //SOC
				{
					if(ipol1==ipol2)
					{
						int spin = ipol1;
						double val = get_onebody_eff_pot(T1, iat1, L1, n1, spin, m1, m2, cal_type, 1);
						VU.at(irc) = complex<double>(val, 0.0);
					}
				}
			}
		}
		vector<complex<double>> potm_tmp(ParaO.nloc, complex<double>(0.0, 0.0));

		// The first term
		pzgemm_(&transN, &transN,
			&NLOCAL, &NLOCAL, &NLOCAL,
			&alpha, 
			VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc,
			VECTOR_TO_PTR(this->Sm_k.at(ik)), &one_int, &one_int, ParaO.desc,
			&beta,
			VECTOR_TO_PTR(potm_tmp), &one_int, &one_int, ParaO.desc);
		
		for(int irc=0; irc<ParaO.nloc; irc++)
		{
			this->pot_eff_k.at(ik).at(irc) += 0.5*potm_tmp.at(irc);
		}


		//The second term
		ZEROS(VECTOR_TO_PTR(potm_tmp), ParaO.nloc);
		
		pzgemm_(&transN, &transN,
			&NLOCAL, &NLOCAL, &NLOCAL,
			&alpha, 
			VECTOR_TO_PTR(this->Sm_k.at(ik)), &one_int, &one_int, ParaO.desc, 
			VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc,
			&beta,
			VECTOR_TO_PTR(potm_tmp), &one_int, &one_int, ParaO.desc);
		
		for(int irc=0; irc<ParaO.nloc; irc++)
		{
			this->pot_eff_k.at(ik).at(irc) += 0.5*potm_tmp.at(irc);
		}
	}
		
	//code for testing whther the effective potential is Hermitian
	/*
	bool pot_Hermitian = true;
	if(GAMMA_ONLY_LOCAL)
	{
		for(int i=0; i<NLOCAL; i++)
		{
			for(int j=0; j<NLOCAL; j++)
			{
				int iic = i*NLOCAL + j;
				int jjc = j*NLOCAL + i;
				double tmp = pot_eff_gamma.at(spin).at(iic) - pot_eff_gamma.at(spin).at(jjc);
				if(tmp>1.0e-9) pot_Hermitian = false;
			}
		}
	}
	else
	{
		for(int i=0; i<NLOCAL; i++)
		{
			for(int j=i; j<NLOCAL; j++)
			{
				int iic = i*NLOCAL + j;
				int jjc = j*NLOCAL + i;
				complex<double> tmp = pot_eff_k.at(ik).at(iic) - conj(pot_eff_k.at(ik).at(jjc));
				double tmp_norm = sqrt(std::norm(tmp));
				if(tmp_norm>1.0e-9) pot_Hermitian = false;
			}
		}
	}


	if(MY_RANK==0)
	{
		ofstream of_potH("Hermitian_pot.dat",ios_base::app);
		of_potH << "Hermitian  " << pot_Hermitian << endl;
	}
	*/
	

	//ofs_running << "dftu.cpp "<< __LINE__  << endl;

	return;	
}


double DFTU::get_onebody_eff_pot
(
	const int T, const int iat,
	const int L, const int N, const int spin, 
	const int m0, const int m1,
    const int type, const bool newlocale 
)
{
	TITLE("DFTU","get_onebody_eff_pot");

	double VU = 0.0;

	//if(!Yukawa && N!=0) return 0.0;

	switch(type)
	{
	case 1:  //rotationally invarient formalism and FLL double counting; Not available at present
		

		break;

	case 2:  //rotationally invarient formalism and AMF double counting; Not available at present
		
		break;

	case 3:	//simplified formalism and FLL double counting
		if(newlocale)
		{
			if(Yukawa)
			{				
				if(m0==m1) VU = (this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*( 0.5 - this->locale.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*this->locale.at(iat).at(L).at(N).at(spin)(m0,m1);
			}
			else
			{
				if(m0==m1) VU = (this->U[T]-this->J[T])*( 0.5 - this->locale.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U[T]-this->J[T])*this->locale.at(iat).at(L).at(N).at(spin)(m0,m1);
			}	
		}
		else
		{	
			if(Yukawa)
			{				
				if(m0==m1) VU = (this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*( 0.5 - this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1);
			}
			else
			{
				if(m0==m1) VU = (this->U[T]-this->J[T])*( 0.5 - this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U[T]-this->J[T])*this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1);
			}
		}
		
		break;

	case 4: //simplified formalism and AMF double counting; ot available at present
		
		break;		
	}

	return VU;
}


void DFTU::force_stress()
{
	TITLE("DFTU", "force_stress");

	if(FORCE)
	{
		for(int iat=0; iat<ucell.nat; iat++)
		{
			for(int dim=0; dim<3; dim++)
			{
				this->force_dftu.at(iat).at(dim) = 0.0;
			}
		}
	}

	if(STRESS)
	{
		for(int dim=0; dim<3; dim++)
		{
			this->stress_dftu.at(dim).at(0) = 0.0;
			this->stress_dftu.at(dim).at(1) = 0.0;
			this->stress_dftu.at(dim).at(2) = 0.0;
		}
	}

    this->allocate_force_stress();

	//=======================================================
	// intialize single particel effective potential matrix
	//=======================================================
	vector<vector<complex<double>>> VU_k;
	vector<vector<double>> VU_gamma;	
	if(FORCE || STRESS)
	{
		this->folding_dSm_soverlap();

		if(GAMMA_ONLY_LOCAL)
		{
			VU_gamma.resize(NSPIN);
			for(int is=0; is<NSPIN; is++)
			{
				VU_gamma.at(is).resize(ParaO.nloc, 0.0);
			}
		}
		else
		{
			if(NSPIN==1)
			{
				VU_k.resize(1);
				VU_k.at(0).resize(ParaO.nloc, complex<double>(0.0,0.0));
			}
			else if(NSPIN==2 || NSPIN==4)
			{
				VU_k.resize(2);
				for(int is=0; is<2; is++)
				{
					VU_k.at(is).resize(ParaO.nloc, complex<double>(0.0,0.0));
				}
			}
		}

		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int T1 = this->iwt2it.at(iwt1);
			const int iat1 = this->iwt2iat.at(iwt1);
			const int L1 = this->iwt2l.at(iwt1);
			const int n1 = this->iwt2n.at(iwt1);
			const int m1 = this->iwt2m.at(iwt1);
			const int ipol1 = this->iwt2ipol.at(iwt1);			

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int T2 = this->iwt2it.at(iwt2);
				const int iat2 = this->iwt2iat.at(iwt2);
				const int L2 = this->iwt2l.at(iwt2);
				const int n2 = this->iwt2n.at(iwt2);
				const int m2 = this->iwt2m.at(iwt2);
				const int ipol2 = this->iwt2ipol.at(iwt2);				

				int irc = ic*ParaO.nrow + ir;

				if(INPUT.orbital_corr[T1]==-1 || INPUT.orbital_corr[T2]==-1) continue;

				if(iat1!=iat2) continue;

				// if(Yukawa)
				// {
					// if(L1<INPUT.orbital_corr[T1] || L2<INPUT.orbital_corr[T2]) continue;
				// }
				// else
				// {
					if(L1!=INPUT.orbital_corr[T1] || L2!=INPUT.orbital_corr[T2] || n1!=0 || n2!=0) continue;
				// }

				if(L1!=L2 || n1!=n2) continue;

				// const int m1_all = m1 + (2*L1+1)*ipol1;
				// const int m2_all = m2 + (2*L2+1)*ipol2;

				if(GAMMA_ONLY_LOCAL)
				{
					for(int is=0; is<NSPIN; is++) //no soc
					{
						double val = get_onebody_eff_pot(T1, iat1, L1, n1, is, m1, m2, cal_type, false);
					
						VU_gamma.at(is).at(irc) = val;
					}
				}
				else
				{
					if(NSPIN==1 || NSPIN==4)
					{
						double val = get_onebody_eff_pot(T1, iat1, L1, n1, 0, m1, m2, cal_type, false);
						VU_k.at(0).at(irc) = complex<double>(val, 0.0);
					}
					else if(NSPIN==2)
					{
						for(int is=0; is<NSPIN; is++)
						{
							double val = get_onebody_eff_pot(T1, iat1, L1, n1, is, m1, m2, cal_type, false);
							VU_k.at(is).at(irc) = complex<double>(val, 0.0);
						}
					}
					else if(NSPIN==4)//SOC
					{
						if(ipol1==ipol2)
						{
							int is = ipol1;
							double val = get_onebody_eff_pot(T1, iat1, L1, n1, is, m1, m2, cal_type, false);
							VU_k.at(is).at(irc) = complex<double>(val, 0.0);
						}
					}
				}					
									
			}//end ic
		}//end ir
	}

	//=======================================================
	//       calculate force
	//=======================================================
	if(FORCE) 
	{
		if(GAMMA_ONLY_LOCAL) cal_force_gamma(VU_gamma);
		else  cal_force_k(VU_k);
	}

	if(STRESS)
	{
		if(GAMMA_ONLY_LOCAL) cal_stress_gamma(VU_gamma);
		else  cal_stress_k(VU_k);
	}			

	return;
}


void DFTU::cal_force_k(vector<vector<complex<double>>> &VU)
{
	TITLE("DFTU", "cal_force_k");

	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;
	
	vector<vector<complex<double>>> ftmp(ucell.nat);
	for(int ia=0; ia<ucell.nat; ia++)
	{
		ftmp.at(ia).resize(3, complex<double>(0.0,0.0)); //three dimension
	}

	vector<vector<complex<double>>> dm_VU_dSm(3);
	for(int dim=0; dim<3; dim++)
	{
		dm_VU_dSm.at(dim).resize(ParaO.nloc, complex<double>(0.0, 0.0));
	}
	
	for(int ik=0; ik<kv.nks; ik++)	
	{
		const int spin = kv.isk[ik];

		for(int dim=0; dim<3; dim++)
		{
			vector<complex<double>> mat_tmp(ParaO.nloc, complex<double>(0.0, 0.0));
			vector<complex<double>> force_tmp(ParaO.nloc, complex<double>(0.0, 0.0));

			if(dim==0) //dim=1,2 are same as dim=0
			{
				ZEROS(VECTOR_TO_PTR(mat_tmp), ParaO.nloc);

				pzgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_k.at(ik).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);
			}			

			pzgemm_(&transN, &transN,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
				this->dSm_k[ik][dim], &one_int, &one_int, ParaO.desc,
				&beta,
				VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) += force_tmp.at(irc);
			}

			//=========================================
			//   the second part
			//=========================================
			ZEROS(VECTOR_TO_PTR(force_tmp), ParaO.nloc);

			pzgemm_(&transN, &transT,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				this->dSm_k[ik][dim], &one_int, &one_int, ParaO.desc, 
				VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
				&beta,
				VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) -= force_tmp.at(irc);
			}
		}//end dim				
	}//end ik

	for(int dim=0; dim<3; dim++)
	{
		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int iat1 = this->iwt2iat.at(iwt1);

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int iat2 = this->iwt2iat.at(iwt2);

				const int irc = ic*ParaO.nrow + ir;

				if(iat1==iat2 && iwt1==iwt2)
				{
					ftmp.at(iat1).at(dim) += dm_VU_dSm.at(dim).at(irc);
				}	
			}//end ic
		}//end ir

		for(int iat=0; iat<ucell.nat; iat++)
		{
			complex<double> val = ftmp.at(iat).at(dim);
			MPI_Allreduce(&val, &ftmp.at(iat).at(dim), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

			this->force_dftu.at(iat).at(dim) = ftmp.at(iat).at(dim).real();
		}
	}

	/*
	if(MY_RANK==0)
	{
		ofstream of_fdftu("force_dftu.dat", ios_base::app);
		for(int iat=0; iat<ucell.nat; iat++)
		{
			of_fdftu << "atom" << iat << "  force_x= " << "(" << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(0).real() << " i" << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(0).imag()
			<< ")  force_y= " << "(" << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(1).real() << " i" << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(1).imag()
			<< ")  force_z= " << "(" << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(2).real() << " i" << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(2).imag() << ")" << endl;								
		}
	}
	*/

	return;
}


void DFTU::cal_stress_k(vector<vector<complex<double>>> &VU)
{
	TITLE("DFTU", "cal_stress_k");

	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;
	
	int count = 0;
	for(int dim1=0; dim1<3; dim1++)
	{
		for(int dim2=dim1; dim2<3; dim2++)
		{
			vector<complex<double>> dm_VU_sover(ParaO.nloc, complex<double>(0.0, 0.0));

			for(int ik=0; ik<kv.nks; ik++)
			{
				const int spin = kv.isk[ik];
				
				// The first term
				vector<complex<double>> stress_tmp(ParaO.nloc, complex<double>(0.0, 0.0));

				//Calculate mat_tmp=dm*VU
				vector<complex<double>> mat_tmp(ParaO.nloc, complex<double>(0.0, 0.0));

				pzgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_k.at(ik).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);

				pzgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					this->soverlap_k[ik][count], &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
				}

				// The second term
				ZEROS(VECTOR_TO_PTR(stress_tmp), ParaO.nloc);

				pzgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					this->soverlap_k[ik][count], &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					dm_VU_sover.at(irc) += 0.5*stress_tmp.at(irc);
				}

			}//end ik

			complex<double> stmp(0.0, 0.0);
			for(int ir=0; ir<ParaO.nrow; ir++)
			{
				const int iwt1 = ParaO.MatrixInfo.row_set[ir];

				for(int ic=0; ic<ParaO.ncol; ic++)
				{
					const int iwt2 = ParaO.MatrixInfo.col_set[ic];

					const int irc = ic*ParaO.nrow + ir;

					if(iwt1==iwt2) stmp += dm_VU_sover.at(irc);

				}//end ic

			}//end ir

			double val = stmp.real();
			MPI_Allreduce(&val, &stress_dftu.at(dim1).at(dim2), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			complex<double> tmp;
			MPI_Allreduce(&stmp, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			/*
			if(MY_RANK==0)
			{
				ofstream of_s("stress_test.dat", ios_base::app);
				
				of_s << "(" << fixed << setprecision(9) << setw(12) << tmp.real() << " i" << fixed << setprecision(9) << setw(12) << tmp.imag() << ")       ";
	
				of_s << endl;
			}
			*/
				
						
			count++;
		}//end dim2
		
	}//end dim1
	
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(i>j) this->stress_dftu.at(i).at(j) = stress_dftu.at(j).at(i);
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			this->stress_dftu.at(i).at(j) *=  ucell.lat0 / ucell.omega;
		}
	}

	/*
	if(MY_RANK==0)
	{
		ofstream of_sdftu("stress_dftu.dat", ios_base::app);
		for(int dim1=0; dim1<3; dim1++)
		{
			for(int dim2=0; dim2<3; dim2++)
			{
				of_sdftu << fixed << setprecision(9) << setw(12) << this->stress_dftu.at(dim1).at(dim2) << "    ";
			}
			of_sdftu << endl;
		}
		of_sdftu << endl;			
	}
	*/

	return;
}


void DFTU::cal_force_gamma(vector<vector<double>> &VU)
{
	TITLE("DFTU", "cal_force_gamma");

	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;
	
	vector<vector<double>> dm_VU_dSm;
	dm_VU_dSm.resize(3);
	for(int dim=0; dim<3; dim++)
	{
		dm_VU_dSm.at(dim).resize(ParaO.nloc, 0.0);
	}

	vector<double> force_tmp(ParaO.nloc, 0.0);
	//Calculate mat_tmp=dm*VU
	vector<double> mat_tmp(ParaO.nloc, 0.0);
	
	for(int ik=0; ik<kv.nks; ik++)
	{
		const int spin = kv.isk[ik];

		for(int dim=0; dim<3; dim++)
		{
			// The first part
			if(dim==0) //dim=1,2 are same as dim=0
			{
				ZEROS(VECTOR_TO_PTR(mat_tmp), ParaO.nloc);

				pdgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_gamma.at(spin).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);
			}

			if(dim==0)
			{
				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					LM.DSloc_x, &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==1)
			{
				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					LM.DSloc_y, &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==2)
			{
				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					LM.DSloc_z, &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			
			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) += force_tmp.at(irc);
			}

			// The second part
			ZEROS(VECTOR_TO_PTR(force_tmp), ParaO.nloc);

			if(dim==0)
			{
				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LM.DSloc_x, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==1)
			{
				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LM.DSloc_y, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==2)
			{
				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LM.DSloc_z, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) -= force_tmp.at(irc);
			}

		}// end dim

	}//end ik

	vector<vector<double>> ftmp;
	ftmp.resize(ucell.nat);
	for(int iat=0; iat<ucell.nat; iat++)
	{
		ftmp.at(iat).resize(3, 0.0);
	}

	for(int dim=0; dim<3; dim++)
	{
		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int iat1 = this->iwt2iat.at(iwt1);

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int iat2 = this->iwt2iat.at(iwt2);

				const int irc = ic*ParaO.nrow + ir;

				if(iat1==iat2 && iwt1==iwt2)
				{
					ftmp.at(iat1).at(dim) += dm_VU_dSm.at(dim).at(irc);
				}	
			}//end ic
		}//end ir
	}//end dim


	for(int iat=0; iat<ucell.nat; iat++)
	{
		for(int dim=0; dim<3; dim++)
		{
			double tmp = 0.0;
			double val = ftmp.at(iat).at(dim);
			MPI_Allreduce(&val, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			this->force_dftu.at(iat).at(dim) = tmp;
		}

		/*
		if(MY_RANK==0)
		{
			ofstream of_fdftu("force_dftu.dat");
			for(int iat=0; iat<ucell.nat; iat++)
			{
				of_fdftu << "atom" << iat << "  force_x= " << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(0)
				<< "  force_y= " << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(1) 
				<< "  force_z= " << fixed << setprecision(9) << setw(12) << ftmp.at(iat).at(2) << endl;								
			}
		}	
		*/
	}

	return;
}


void DFTU::cal_stress_gamma(vector<vector<double>> &VU)
{
	TITLE("DFTU", "cal_stress_gamma");

	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;
	
	int count = 0;
	for(int dim1=0; dim1<3; dim1++)
	{
		for(int dim2=dim1; dim2<3; dim2++)
		{
			vector<double> dm_VU_sover(ParaO.nloc, 0.0);

			for(int ik=0; ik<kv.nks; ik++)
			{
				const int spin = kv.isk[ik];
				
				// The first term
				vector<double> stress_tmp(ParaO.nloc, 0.0);

				//Calculate mat_tmp=dm*VU
				vector<double> mat_tmp(ParaO.nloc, 0.0);

				pdgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_gamma.at(spin).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);

				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					this->soverlap_gamma[count], &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
				}

				// The second term
				ZEROS(VECTOR_TO_PTR(stress_tmp), ParaO.nloc);

				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					this->soverlap_gamma[count], &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					dm_VU_sover.at(irc) += 0.5*stress_tmp.at(irc);
				}

			}//end ik

			double stmp = 0.0;
			for(int ir=0; ir<ParaO.nrow; ir++)
			{
				const int iwt1 = ParaO.MatrixInfo.row_set[ir];

				for(int ic=0; ic<ParaO.ncol; ic++)
				{
					const int iwt2 = ParaO.MatrixInfo.col_set[ic];

					const int irc = ic*ParaO.nrow + ir;

					if(iwt1==iwt2) stmp += dm_VU_sover.at(irc);

				}//end ic

			}//end ir

			MPI_Allreduce(&stmp, &stress_dftu.at(dim1).at(dim2), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
			
			count++;
		}//end dim2
		
	}//end dim1

	
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(i>j) this->stress_dftu.at(i).at(j) = stress_dftu.at(j).at(i);
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			this->stress_dftu.at(i).at(j) *=  ucell.lat0 / ucell.omega;
		}
	}

	/*
	if(MY_RANK==0)
	{
		ofstream of_fdftu("force_dftu.dat");
		for(int iat=0; iat<ucell.nat; iat++)
		{
			of_fdftu << "atom" << iat << "  force_x=" << "" << fixed << setprecision(9) << setw(12) << this->force_dftu.at(iat).at(0)
			<< "  force_y=" <<  fixed << setprecision(9) << setw(12) << this->force_dftu.at(iat).at(1) 
			<< "  force_z=" << fixed << setprecision(9) << setw(12) << this->force_dftu.at(iat).at(2) << endl;								
		}

		if(STRESS)
		{
			ofstream of_sdftu("stress_dftu.dat");
			for(int dim0=0; dim0<3; dim0++)
			{
				for(int dim1=0; dim1<3; dim1++)
				{
					of_sdftu << fixed << setprecision(9) << setw(12) << this->stress_dftu.at(dim0).at(dim1) << "     ";
				}
				of_sdftu << endl;
			}
		}
	}
	*/

	return;
}


void DFTU::folding_dSm_soverlap()
{
	TITLE("DFTU", "folding_dSm_soverlap");

	int nnr = 0;

	if(GAMMA_ONLY_LOCAL)
	{
		if(STRESS)
		{
			for(int i=0; i<6; i++)
			{
				ZEROS(soverlap_gamma[i], ParaO.nloc);
			}			
		}
	}
	else
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			for(int dim=0; dim<3; dim++)
			{
				ZEROS(dSm_k[ik][dim], ParaO.nloc);
			}
		}

		if(STRESS)
		{
			for(int ik=0; ik<kv.nks; ik++)
			{
				for(int i=0; i<6; i++)
				{
					ZEROS(soverlap_k[ik][i], ParaO.nloc);
				}
			}
		}
	}
	

	Vector3<double> tau1, tau2, dtau;
	Vector3<double> dtau1, dtau2, tau0;
    for(int T1=0; T1<ucell.ntype; ++T1)
    {
		Atom* atom1 = &ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; ++I1)
        {
			tau1 = atom1->tau[I1];
            
            GridD.Find_atom(tau1, T1, I1);
            for(int ad=0; ad<GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);

				Atom* atom2 = &ucell.atoms[T2];

				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;

				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();				

				if(distance < rcut)
				{
					int iw1_all = ucell.itiaiw2iwt( T1, I1, 0) ; //iw1_all = combined index (it, ia, iw)

					for(int jj=0; jj<atom1->nw*NPOL; ++jj)
					{
						const int jj0 = jj/NPOL;
						const int L1 = atom1->iw2l[jj0];
						const int N1 = atom1->iw2n[jj0];
						const int m1 = atom1->iw2m[jj0];
						int iw2_all = ucell.itiaiw2iwt( T2, I2, 0);

						for(int kk=0; kk<atom2->nw*NPOL; ++kk)
						{
							const int kk0 = kk/NPOL;
							const int L2 = atom2->iw2l[kk0];
							const int N2 = atom2->iw2n[kk0];
							const int m2 = atom2->iw2m[kk0];
							
							if ( !ParaO.in_this_processor(iw1_all,iw2_all) )
							{
								++iw2_all;
								continue;
							}

							int mu = ParaO.trace_loc_row[iw1_all];
							int nu = ParaO.trace_loc_col[iw2_all];
							int irc = nu*ParaO.nrow + mu;
														
							if(GAMMA_ONLY_LOCAL)
							{
								if(STRESS)
								{
									this->soverlap_gamma[0][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+0];
									this->soverlap_gamma[1][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+1];
									this->soverlap_gamma[2][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+2];
									this->soverlap_gamma[3][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+1];
									this->soverlap_gamma[4][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+2];
									this->soverlap_gamma[5][irc] += LM.DSloc_Rz[nnr]*LM.DH_r[nnr*3+2];
								}
							}
							else
							{
								Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z); 
							
								for(int ik=0; ik<kv.nks; ik++)
								{								
									const double arg = ( kv.kvec_d[ik] * dR ) * TWO_PI;
									const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );

									this->dSm_k[ik][0][irc] += LM.DSloc_Rx[nnr]*kphase;
									this->dSm_k[ik][1][irc] += LM.DSloc_Ry[nnr]*kphase;
									this->dSm_k[ik][2][irc] += LM.DSloc_Rz[nnr]*kphase;

									if(STRESS)
									{																												
										this->soverlap_k[ik][0][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+0]*kphase;
										this->soverlap_k[ik][1][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+1]*kphase;
										this->soverlap_k[ik][2][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+2]*kphase;
										this->soverlap_k[ik][3][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+1]*kphase;
										this->soverlap_k[ik][4][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+2]*kphase;
										this->soverlap_k[ik][5][irc] += LM.DSloc_Rz[nnr]*LM.DH_r[nnr*3+2]*kphase;																
									}
								}	
							}
																																																																				
							++nnr;													
							++iw2_all;
						}// nw2 

						++iw1_all;
						
					}// nw1
				}// distance
				else if(distance>=rcut)
				{
					int start1 = ucell.itiaiw2iwt( T1, I1, 0);
					int start2 = ucell.itiaiw2iwt( T2, I2, 0);
					bool is_adj = false;
					for (int ad0=0; ad0<GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						
						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();
						if(distance1<rcut1 && distance2<rcut2)
						{
							is_adj = true;
							break;
						}
					}//ad0
					if( is_adj )
					{
						for(int jj=0; jj<atom1->nw * NPOL; ++jj)
						{
							const int mu = ParaO.trace_loc_row[start1+jj];
							if(mu<0)continue; 

							for(int kk=0; kk<atom2->nw * NPOL; ++kk)
							{
								const int nu = ParaO.trace_loc_col[start2+kk];
								if(nu<0)continue;

								++nnr;
							}//kk
						}//jj
					}
				}//distance
			}// ad
		}// I1
	}// T1

	return;
}


void DFTU::allocate_force_stress()
{
	TITLE("DFTU","allocate_force_stress");

	if(GAMMA_ONLY_LOCAL)
	{
		if(STRESS)
		{
			// this->soverlap_gamma.resize(6);   //xx, xy, xz, yy, yz, zz
			// for(int i=0; i<6; i++)
			// {				
				// this->soverlap_gamma.at(i).resize(ParaO.nloc, 0.0);
			// }
			soverlap_gamma = new double* [6];
			for(int i=0; i<6; i++)
			{				
				soverlap_gamma[i] = new double [ParaO.nloc];
			}
		}
	}
	else
	{
		dSm_k = new complex<double>** [kv.nks];
		//this->dSm_k.resize(kv.nks);

		for(int ik=0; ik<kv.nks; ik++)
		{
			//this->dSm_k.at(ik).resize(3);
			dSm_k[ik] = new complex<double>* [3];
			for(int dim=0; dim<3; dim++)
			{
				//this->dSm_k.at(ik).at(dim).resize(ParaO.nloc, ZERO);
				dSm_k[ik][dim] = new complex<double>[ParaO.nloc];
			}
		}

		if(STRESS)
		{
			//this->soverlap_k.resize(kv.nks);
			soverlap_k = new complex<double>** [kv.nks];

			for(int ik=0; ik<kv.nks; ik++)
			{
				//this->soverlap_k.at(ik).resize(6);   //xx, xy, xz, yy, yz, zz
				soverlap_k[ik] = new complex<double>* [6];
				for(int i=0; i<6; i++)
				{				
					//this->soverlap_k.at(ik).at(i).resize(ParaO.nloc, ZERO);
					soverlap_k[ik][i] = new complex<double> [ParaO.nloc];
				}
			}
		}
	}
	
	return;
}


void DFTU::erase_force_stress()
{
	TITLE("DFTU","erase_force_stress");

	if(GAMMA_ONLY_LOCAL)
	{
		if(STRESS)
		{
			//this->soverlap_gamma.resize(6);   //xx, xy, xz, yy, yz, zz
			for(int i=0; i<6; i++)
			{				
				//this->soverlap_gamma.at(i).resize(1, 0.0);
				delete [] soverlap_gamma[i];
			}
			delete [] soverlap_gamma;

		}
	}
	else
	{
		//this->dSm_k.resize(kv.nks);

		for(int ik=0; ik<kv.nks; ik++)
		{
			//this->dSm_k.at(ik).resize(3);
			for(int dim=0; dim<3; dim++)
			{
				//this->dSm_k.at(ik).at(dim).resize(1, ZERO);
				delete [] dSm_k[ik][dim];
			}
			delete [] dSm_k[ik];
		}
		delete [] dSm_k;

		if(STRESS)
		{
			//this->soverlap_k.resize(kv.nks);

			for(int ik=0; ik<kv.nks; ik++)
			{
				//this->soverlap_k.at(ik).resize(6);   //xx, xy, xz, yy, yz, zz
				for(int i=0; i<6; i++)
				{		
					//this->soverlap_k.at(ik).at(i).resize(1, ZERO);
					delete [] soverlap_k[ik][i];
				}
				delete [] soverlap_k[ik];
			}
			delete [] soverlap_k;
		}
	}
			
	return;
}


void DFTU::cal_yukawa_lambda()
{
	TITLE("DFTU", "cal_yukawa_lambda");
	
	double sum_rho = 0.0;
	double sum_rho_lambda = 0.0;	
	for(int is=0; is<NSPIN; is++)
	{
		for(int ir=0; ir<pw.nrxx; ir++) 
		{
			double rho_ir = chr.rho[is][ir];
			sum_rho += rho_ir;

			double lambda_ir = 2*pow(3*rho_ir/PI, (double)1.0/6.0);
			sum_rho_lambda += lambda_ir*rho_ir;
		}
	}

	double val1 = 0.0;
	double val2 = 0.0;
	MPI_Allreduce(&sum_rho, &val1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_rho_lambda, &val2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	this->lambda = val2/val1;

	double lambda1 = this->lambda;


	//rescaling
	this->lambda /= 1.6;

	/*
	double Nval = 0.0;
	double ele_cor = 0.0;
	for(int T=0; T<ucell.ntype; T++)
	{
		const int L = INPUT.orbital_corr[T];
		const int N = ucell.atoms[T].l_nchi[L];
		Nval += ucell.atoms[T].na*ucell.atoms[T].zv;

		if(L==-1) continue;	

		for(int I=0; I<ucell.atoms[T].na; I++)
		{
			const int iat = ucell.itia2iat(T, I);
			
			for(int n=0; n<N; n++)
			{
				for(int m=0; m<2*L+1; m++)
				{					
					for(int is=0; is<2; is++)
					{
						ele_cor += this->locale.at(iat).at(L).at(n).at(is)(m,m);
					}											
				}
			}
		}
	}
	*/

	double rho_remain = (Nval-Nc)/ucell.omega;
	double lambda2 = 2*pow(3*rho_remain/PI, (double)1.0/6.0);
	
	/*
	if(MY_RANK==0)
	{
		ofstream ofs_lambda("lambda.dat", ios_base::app);

		ofs_lambda << "All valence electrons.   lambda=" << fixed << setw(8) << setprecision(4) << lambda1 << endl;
		ofs_lambda << "Subtracting correltaed electrons from valence electrons.   lambda=" << fixed << setw(8) << setprecision(4) << lambda2 << endl;
		ofs_lambda << endl;
	}
	*/


	return;
}


void DFTU::cal_slater_Fk(const int L, const int T)
{
	TITLE("DFTU","cal_slater_Fk");
	
	if(Yukawa)
	{	 
		//this->lambda = INPUT.yukawa_lambda;
		
		for(int chi=0; chi<ucell.atoms[T].l_nchi[L]; chi++)
		{
		 //	if(chi!=0) continue;
			const int mesh = ORB.Phi[T].PhiLN(L,chi).getNr();
	
			for(int k=0; k<=L; k++)
			{			
				for(int ir0=1; ir0<mesh; ir0++)
				{
					double r0 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir0);
					const double rab0 = ORB.Phi[T].PhiLN(L,chi).getRab(ir0);
					const double R_L0 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir0);
	
					for(int ir1=1; ir1<mesh; ir1++) 
					{
						double bslval, hnkval;
						double r1 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir1);
						const double rab1 = ORB.Phi[T].PhiLN(L,chi).getRab(ir1);
						const double R_L1 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir1);		
						
						int l = 2*k;
						if(ir0<ir1)  //less than
						{
						 	sphbsl_(&l, &r0, &lambda, &bslval);
							sphhnk_(&l, &r1, &lambda, &hnkval);
						}
						else //greater than
						{
						 	sphbsl_(&l, &r1, &lambda, &bslval);
							sphhnk_(&l, &r0, &lambda, &hnkval);
						}					
						this->Fk.at(T).at(L).at(chi).at(k) -= (4*k+1)*lambda*pow(R_L0,2)*bslval*hnkval*pow(R_L1,2)*pow(r0,2)*pow(r1,2)*rab0*rab1;					
					}
				}
			 //	this->Fk.at(T).at(chi).at(k) /= pow(norm,2);
			}			
		}		
	}

	return;
}


void DFTU::cal_slater_UJ(const int istep, const int iter)
{
	TITLE("DFTU", "cal_slater_UJ");
	if(!Yukawa) return;

	this->cal_yukawa_lambda();
	
	for(int it=0; it<ucell.ntype; it++)
	{			
		const int NL = ucell.atoms[it].nwl + 1;
	
		for(int l=0; l<NL; l++)
		{
			int N = ucell.atoms[it].l_nchi[l];
			for(int n=0; n<N; n++)
			{
				ZEROS(VECTOR_TO_PTR(this->Fk.at(it).at(l).at(n)), l+1);
			}			
		}			 	
	}
 	
	for(int T=0; T<ucell.ntype; T++)
	{			
		const int NL = ucell.atoms[T].nwl + 1;

		for(int L=0; L<NL; L++)
		{
			const int N = ucell.atoms[T].l_nchi[L];

			if(L>=INPUT.orbital_corr[T] && INPUT.orbital_corr[T]!=-1)
			{
 				if(L!=INPUT.orbital_corr[T]) continue;

				this->cal_slater_Fk(L, T);

				for(int n=0; n<N; n++)
				{
 					if(n!=0) continue;

					switch(L)
					{
					case 1: //p electrons
    				    if(Yukawa)
						{
							this->U_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(0);
    				    	this->J_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(1)/5.0;
						}
						else
						{
						 //	if(n!=0) continue;

							this->Fk.at(T).at(L).at(n).at(0) = this->U[T];
    				    	this->Fk.at(T).at(L).at(n).at(1) = 5.0*this->J[T];
						}			
						break;

    				case 2: //d electrons
    				    if(Yukawa)
						{
							this->U_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(0);
    				    	this->J_Yukawa.at(T).at(L).at(n) = (this->Fk.at(T).at(L).at(n).at(1)+this->Fk.at(T).at(L).at(n).at(2))/14.0;
						}
						else
						{
						 //	if(n!=0) continue;

							this->Fk.at(T).at(L).at(n).at(0) = this->U[T];
							this->Fk.at(T).at(L).at(n).at(1) = 14.0*this->J[T]/(1.0+0.625);
							this->Fk.at(T).at(L).at(n).at(2) = 0.625*this->Fk.at(T).at(L).at(n).at(1);
						}
						break;

    				case 3: //f electrons
    				    if(Yukawa)
						{
							this->U_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(0);
    				    	this->J_Yukawa.at(T).at(L).at(n) = (286.0*this->Fk.at(T).at(L).at(n).at(1) +
									195.0*this->Fk.at(T).at(L).at(n).at(2)+250.0*this->Fk.at(T).at(L).at(n).at(3))/6435.0;
						}
						else
						{
						 //	if(n!=0) continue;

							this->Fk.at(T).at(L).at(n).at(0) = this->U[T];
							this->Fk.at(T).at(L).at(n).at(1) = 6435.0*this->J[T]/(286.0+195.0*0.668+250.0*0.494);
							this->Fk.at(T).at(L).at(n).at(2) = 0.668*this->Fk.at(T).at(L).at(n).at(1);
							this->Fk.at(T).at(L).at(n).at(3) = 0.494*this->Fk.at(T).at(L).at(n).at(1);
						}
						break;
					}

					//Hartree to Rydeberg
					this->U_Yukawa.at(T).at(L).at(n) *= 2.0;
					this->J_Yukawa.at(T).at(L).at(n) *= 2.0;
				}//end n

			}//end if

		}//end L			 	
	}// end T

	/*
	if(MY_RANK==0 && Yukawa)
	{
		ofstream of_UJ("Yukawa_UJ.dat", ios_base::app);
		of_UJ << "ISTEP= " << istep << "  ITER= " << iter << endl;
		of_UJ << "Lambda= " << this->lambda << endl;

		for(int T=0; T<ucell.ntype; T++)
		{			
			const int NL = ucell.atoms[T].nwl + 1;

			for(int L=0; L<NL; L++)
			{
				const int N = ucell.atoms[T].l_nchi[L];

				if(L>=INPUT.orbital_corr[T] && INPUT.orbital_corr[T]!=-1)
				{
 if(L!=INPUT.orbital_corr[T]) continue;
					for(int n=0; n<N; n++)
					{
 if(n!=0) continue;						
						double Ueff = (this->U_Yukawa.at(T).at(L).at(n) - this->J_Yukawa.at(T).at(L).at(n))*Ry_to_eV;
						of_UJ << "atom type=" << T << "  L=" << L << "  chi=" << n << endl;
						of_UJ << "U[" << n << "]=" << this->U_Yukawa.at(T).at(L).at(n)*Ry_to_eV << "    " << "J[" << n << "]=" << this->J_Yukawa.at(T).at(L).at(n)*Ry_to_eV
						<< "    Ueff[" << n << "]=" << Ueff << endl;
						of_UJ << endl; 
					}
				}
			}
		}
	}
	*/
		
	return;
}


void DFTU::output()
{
	TITLE("DFTU", "output");
	
	ofs_running << "//=========================L(S)DA+U===========================//" << endl;

	for(int T=0; T<ucell.ntype; T++)
	{			
		const int NL = ucell.atoms[T].nwl + 1;

		for(int L=0; L<NL; L++)
		{
			const int N = ucell.atoms[T].l_nchi[L];

			if(L>=INPUT.orbital_corr[T] && INPUT.orbital_corr[T]!=-1)
			{
 				if(L!=INPUT.orbital_corr[T]) continue;
				
				if(!Yukawa)
				{
					ofs_running << "atom_type=" << T << "  L=" << L << "  chi=" << 0 << "    U=" << this->U[T]*Ry_to_eV 
					<< "ev    J=" << this->J[T]*Ry_to_eV << "ev" << endl;
				}
				else
				{
					for(int n=0; n<N; n++)
					{
 						if(n!=0) continue;						
						double Ueff = (this->U_Yukawa.at(T).at(L).at(n) - this->J_Yukawa.at(T).at(L).at(n))*Ry_to_eV;
						ofs_running << "atom_type=" << T << "  L=" << L << "  chi=" << n << "    U=" << this->U_Yukawa.at(T).at(L).at(n)*Ry_to_eV << "ev    " << "J=" << this->J_Yukawa.at(T).at(L).at(n)*Ry_to_eV
						<< "ev" << endl;
					}
				}
			}
		}
	}

	ofs_running << "Local occupation matrices" << endl;
	for(int T=0; T<ucell.ntype; T++)
	{
		if(INPUT.orbital_corr[T]==-1) continue;
		const int NL = ucell.atoms[T].nwl+1;
		const int LC = INPUT.orbital_corr[T];

		for(int I=0; I<ucell.atoms[T].na; I++)
		{
			const int iat = ucell.itia2iat(T, I);
			ofs_running << "atoms" << " " << iat ;

			for(int l=0; l<NL; l++)
			{
				// if(Yukawa)
				// {
					// if(l<INPUT.orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=INPUT.orbital_corr[T]) continue;
				// }
				if(l!=INPUT.orbital_corr[T]) continue;
				
				const int N = ucell.atoms[T].l_nchi[l];
				ofs_running << "   L" << " " << l ;

				for(int n=0; n<N; n++)		
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

					ofs_running << "   zeta" << " " << n << endl;
					
					for(int is=0; is<2; is++)
					{
						ofs_running << "spin" << "  " << is << endl;
						for(int m0=0; m0<2*l+1; m0++)
						{
							for(int m1=0; m1<2*l+1; m1++)
							{
								ofs_running << fixed << setw(12) << setprecision(8) << locale.at(iat).at(l).at(n).at(is)(m0, m1);
							}
							ofs_running << endl;
						}
					}					
				}//n
			}//l
		}//I
	}//T
	ofs_running << "//=======================================================//" << endl;
	
	return;
}

/*
void DFTU::print(const int T, const int iat, const int L, const int N, const int iter)
{
	if(MY_RANK) return;
	stringstream st;
	st << "loc_occup_m" << iat << ".dat";
	ofstream ofs_dftu(st.str(), ios_base::app);

	ofs_dftu << "Iter  " << iter << endl;

	const int M = 2*L + 1;
 //int n = N;
	for(int n=0; n<N; n++)
	{
		if(!Yukawa && n!=0) continue;
		// if(n!=0) continue;
		ofs_dftu << "zeta="	<< n << endl;
	 //	ofs_pot << "zeta="	<< n << endl;
	 	double tot_elec = 0.0;

		for(int spin=0; spin<2; spin++)
		{
			ofs_dftu << "Spin  " << spin << endl;
	//		ofs_pot << "Spin  " << spin << endl;
			for(int m1=0; m1<M; m1++)
			{
				for(int m2=0; m2<M; m2++)
				{
					ofs_dftu << fixed << setw(12) << setprecision(8) << locale.at(iat).at(L).at(n).at(spin)(m1,m2);
					if(m1==m2) tot_elec += locale[iat][L][n][spin](m1,m2);

				//	double val = get_onebody_eff_pot(T, iat, L, n, spin, m1, m2, cal_type, 1);
				//	val *= Hartree_to_eV;
				//	ofs_pot << fixed << setw(12) << setprecision(9) << val << "   ";

				}
				ofs_dftu << endl;
			//	ofs_pot << endl;
			}
		}
    // this->Nc = tot_elec;
		ofs_dftu << "Total number of electrons (" << L << "," << n << ")    " << tot_elec << endl;
	}

	ofs_dftu << endl;
//	ofs_pot << endl;

	return;
}


void DFTU::cal_unscreened_slater_Fk(const int L, const int T)
{
	TITLE("DFTU","cal_slater_Fk");

	for(int chi=0; chi<ucell.atoms[T].l_nchi[L]; chi++)
	{
		const int mesh = ORB.Phi[T].PhiLN(L,chi).getNr();

		for(int k=0; k<=L; k++)
		{			
			for(int ir0=1; ir0<mesh; ir0++)
			{
				double r0 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir0);
				const double rab0 = ORB.Phi[T].PhiLN(L,chi).getRab(ir0);
				const double R_L0 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir0);

				for(int ir1=1; ir1<mesh; ir1++) 
				{
					double numerator, denominator;
					double r1 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir1);
					const double rab1 = ORB.Phi[T].PhiLN(L,chi).getRab(ir1);
					const double R_L1 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir1);		
					
					int l = 2*k;
					if(ir0<ir1)  //less than
					{
					 	numerator  = pow(r0, l);
						denominator = pow(r1, l+1);
					}
					else //greater than
					{
					 	numerator  = pow(r1, l);
						denominator = pow(r0, l+1);
					}					
					this->Fk.at(T).at(chi).at(k) += (numerator/denominator)*pow(R_L0,2)*pow(R_L1,2)*pow(r0,2)*pow(r1,2)*rab0*rab1;					
				}
			}
		}
		
		
	}

	this->cal_slater_UJ(T, L);

	return;
}

void DFTU::cal_slater_Vsc(const int T, const int L)
{
	TITLE("DFTU", "cal_slater_Vsc");

	for(int N=0; N<ucell.atoms[T].l_nchi[L]; N++)
	{
		if(!Yukawa && N!=0) continue;

		for(int m0=0; m0<2*L+1; m0++)
		{
			const int gindex0 = L*L + m0;

			for(int m1=0; m1<2*L+1; m1++)
			{				
				const int gindex1 = L*L + m1;

				for(int m2=0; m2<2*L+1; m2++)
				{
					const int gindex2 = L*L + m2;
					const int M0 = m0*(2*L+1) + m2;

					for(int m3=0; m3<2*L+1; m3++)
					{
						const int M1 = m1*(2*L+1) + m3;
						const int gindex3 = L*L + m3;

						for(int k=0; k<=2*L; k+=2)
						{
							const int l = (int)(k/2);
							for(int q=0; q<2*k+1; q++)
							{
								int gindex = k*k + q;

								double gaunt1 = UOT.get_Gaunt_coefficients(gindex0, gindex2, gindex);
								double gaunt2 = UOT.get_Gaunt_coefficients(gindex1, gindex3, gindex);

								this->Vsc.at(T).at(N)(M0, M1) += FOUR_PI*gaunt1*gaunt2*Fk.at(T).at(N).at(l)/(2.0*k+1.0);
							}
						}
					}
				}
			}
		}
	}

	return;
}
*/
