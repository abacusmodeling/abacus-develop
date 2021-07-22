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
#include "../module_base/constants.h"
#include "../src_pw/global.h"
#include "global_fp.h"
#include "../module_base/global_function.h"
#include "../module_base/inverse_matrix.h"
#include "LOOP_ions.h"
#include "LCAO_matrix.h"
#include "../src_pw/magnetism.h"
#include "../module_orbital/ORB_gen_tables.h"
#include "../src_pw/charge.h"
#include "LCAO_nnr.h"

extern "C"
{
  void pzgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const std::complex<double> *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const std::complex<double> *beta,
		std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
  
  void pdgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC);

  void pztranc_(
    const int *M, const int *N,
    const std::complex<double> *alpha,
    const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
    const std::complex<double> *beta,
    std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);

  void pdtran_(
    const int *M, const int *N,
    const double *alpha,
    const double *A, const int *IA, const int *JA, const int *DESCA,
    const double *beta,
    double *C, const int *IC, const int *JC, const int *DESCC);
}

DFTU dftu;

DFTU::DFTU(){}

DFTU::~DFTU(){}

void DFTU::init(
	UnitCell_pseudo &cell, // unitcell class
	Parallel_Orbitals &po // parallel orbitals parameters
)
{
    TITLE("DFTU", "init");

	// please, do not use 'INPUT' directly in the class!
	// needs reconstructions in future
	// global parameters, need to be removed in future
	const int npol = NPOL; // number of polarization directions
	const int nlocal = NLOCAL; // number of total local orbitals
	const int nks = kv.nks; // number of k-points
	const int nspin = NSPIN; // number of spins
	const int dftu_type = INPUT.dftu_type;
	const int double_counting = INPUT.double_counting;
	

	if(dftu_type==1 && double_counting==1) cal_type = 1;
	else if(dftu_type==1 && double_counting==2) cal_type = 2;
	else if(dftu_type==2 && double_counting==1) cal_type = 3;
	else if(dftu_type==2 && double_counting==2) cal_type = 4;
	else WARNING_QUIT("DFT+U", "Wrong parameter");
		
	this->EU = 0.0;

	if(FORCE)
	{
		this->force_dftu.resize(cell.nat);
		for(int ia=0; ia<cell.nat; ia++)
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

	this->locale.resize(cell.nat);
	this->locale_save.resize(cell.nat);

	this->iatlnmipol2iwt.resize(cell.nat);
	this->iat2it.resize(cell.nat);
	this->iwt2it.resize(nlocal);
	this->iwt2iat.resize(nlocal);
	this->iwt2l.resize(nlocal);
	this->iwt2n.resize(nlocal);
	this->iwt2m.resize(nlocal);
	this->iwt2ipol.resize(nlocal);	

	for(int i=0; i<nlocal; i++)
	{
		this->iwt2it.at(i) = -1;
		this->iwt2iat.at(i) = -1;
		this->iwt2l.at(i) = -1;
		this->iwt2n.at(i) = -1;
		this->iwt2m.at(i) = -1;
		this->iwt2ipol.at(i) = -1;
	}

	for(int it=0; it<cell.ntype; ++it)
	{				
		for(int ia=0; ia<cell.atoms[it].na; ia++)
		{
			const int iat = cell.itia2iat(it, ia);
			this->iat2it.at(iat) = it;

			locale.at(iat).resize(cell.atoms[it].nwl+1);
			locale_save.at(iat).resize(cell.atoms[it].nwl+1);

			for(int l=0; l<=cell.atoms[it].nwl; l++)
			{			
				const int N = cell.atoms[it].l_nchi[l];

				locale.at(iat).at(l).resize(N);
				locale_save.at(iat).at(l).resize(N);

				for(int n=0; n<N; n++)
				{
					if(nspin==1 || nspin==2)
					{
						locale.at(iat).at(l).at(n).resize(2);
						locale_save.at(iat).at(l).at(n).resize(2);

						locale.at(iat).at(l).at(n).at(0).create(2*l+1, 2*l+1);
						locale.at(iat).at(l).at(n).at(1).create(2*l+1, 2*l+1);

						locale_save.at(iat).at(l).at(n).at(0).create(2*l+1, 2*l+1);
						locale_save.at(iat).at(l).at(n).at(1).create(2*l+1, 2*l+1);
					}
					else if(nspin==4) //SOC
					{
						locale.at(iat).at(l).at(n).resize(1);
						locale_save.at(iat).at(l).at(n).resize(1);

						locale.at(iat).at(l).at(n).at(0).create((2*l+1)*npol, (2*l+1)*npol);
						locale_save.at(iat).at(l).at(n).at(0).create((2*l+1)*npol, (2*l+1)*npol);
					}												
				}
			}

			//initialize the arrry iatlnm2iwt[iat][l][n][m]
			this->iatlnmipol2iwt.at(iat).resize(cell.atoms[it].nwl+1);
			for(int L=0; L<=cell.atoms[it].nwl; L++)
			{
				this->iatlnmipol2iwt.at(iat).at(L).resize(cell.atoms[it].l_nchi[L]);

				for(int n=0; n<cell.atoms[it].l_nchi[L]; n++)
				{
					this->iatlnmipol2iwt.at(iat).at(L).at(n).resize(2*L+1);
					
					for(int m=0; m<2*L+1; m++)
					{
						this->iatlnmipol2iwt.at(iat).at(L).at(n).at(m).resize(npol);
					}
				}
			}

			for(int iw=0; iw<cell.atoms[it].nw*npol; iw++)
			{
				int iw0 = iw/npol;
				int ipol = iw%npol;
				int iwt = cell.itiaiw2iwt(it, ia, iw);
				int l = cell.atoms[it].iw2l[iw0];
				int n = cell.atoms[it].iw2n[iw0];
				int m = cell.atoms[it].iw2m[iw0];
								
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
		this->Fk.resize(cell.ntype);
		
		this->U_Yukawa.resize(cell.ntype);
		this->J_Yukawa.resize(cell.ntype);	

		for(int it=0; it<cell.ntype; it++)
		{			
			const int NL = cell.atoms[it].nwl + 1;

			this->Fk.at(it).resize(NL);		
			this->U_Yukawa.at(it).resize(NL);
			this->J_Yukawa.at(it).resize(NL);	

			for(int l=0; l<NL; l++)
			{
				int N = cell.atoms[it].l_nchi[l];

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
					if(NSPIN==4)
					{
						locale_save[iat][l][n][0] = locale[iat][l][n][0];

            locale[iat][l][n][0].zero_out();
					}
					else if(NSPIN==1 || NSPIN==2)
					{
						locale_save[iat][l][n][0] = locale[iat][l][n][0];
						locale_save[iat][l][n][1] = locale[iat][l][n][1];

            locale[iat][l][n][0].zero_out();
            locale[iat][l][n][1].zero_out();
					}
				}
			}			
		}
	}

	//=================Part 1======================
	//call PBLAS routine to calculate the product of the S and density matrix
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const complex<double> alpha(1.0,0.0), beta(0.0,0.0);

	vector<complex<double>> srho(ParaO.nloc);
  vector<complex<double>> Sk(ParaO.nloc);
	
	for(int ik=0; ik<kv.nks; ik++)
	{
		// srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
    this->folding_overlap_matrix(ik, &Sk[0]);

		pzgemm_(&transN, &transT,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				&Sk[0], &one_int, &one_int, ParaO.desc, 
				LOC.wfc_dm_2d.dm_k.at(ik).c, &one_int, &one_int, ParaO.desc,
				&beta, 
				&srho[0], &one_int, &one_int, ParaO.desc);

    const int spin = kv.isk[ik];
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

		  			//Calculate the local occupation number matrix			
		  			for(int m0=0; m0<2*l+1; m0++)
		  			{
		  				for(int ipol0=0; ipol0<NPOL; ipol0++)
		  				{
		  					const int iwt0 = this->iatlnmipol2iwt[iat][l][n][m0][ipol0];
		  					const int mu = ParaO.trace_loc_row[iwt0];
		  					const int mu_prime = ParaO.trace_loc_col[iwt0];

		  					for(int m1=0; m1<2*l+1; m1++)
		  					{
		  						for(int ipol1=0; ipol1<NPOL; ipol1++)
		  						{									
		  							const int iwt1 = this->iatlnmipol2iwt[iat][l][n][m1][ipol1];
		  							const int nu = ParaO.trace_loc_col[iwt1];
		  							const int nu_prime = ParaO.trace_loc_row[iwt1];

		  							const int irc = nu*ParaO.nrow + mu;
		  							const int irc_prime = mu_prime*ParaO.nrow + nu_prime;

		  							const int m0_all = m0 + ipol0*(2*l+1);
		  							const int m1_all = m1 + ipol1*(2*l+1);

		  							if( (nu>=0) && (mu>=0) )
		  								locale[iat][l][n][spin](m0_all, m1_all) += (srho[irc]).real()/4.0;									

		  							if( (nu_prime>=0) && (mu_prime>=0) )
		  								locale[iat][l][n][spin](m0_all, m1_all) += (std::conj(srho[irc_prime])).real()/4.0;
		  						}//ipol1										
		  					}//m1
		  				}//ipol0
		  			}//m0
		  		}//end n
		  		// this->print(it, iat, l, N, iter);
		  	}//end l
		  }//end ia
	  }//end it
	}//ik

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

					if(NSPIN==1 || NSPIN==4)
					{
            matrix temp(locale[iat][l][n][0]);
						MPI_Allreduce( &temp(0,0), &locale[iat][l][n][0](0,0), (2*l+1)*NPOL*(2*l+1)*NPOL,
										      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
					else if(NSPIN==2)
					{
            matrix temp0(locale[iat][l][n][0]);
						MPI_Allreduce( &temp0(0,0), &locale[iat][l][n][0](0,0), (2*l+1)*(2*l+1),
										      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

						matrix temp1(locale[iat][l][n][1]);
						MPI_Allreduce( &temp1(0,0), &locale[iat][l][n][1](0,0), (2*l+1)*(2*l+1),
										      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
				
					// for the case spin independent calculation
					switch(NSPIN)
					{
					  case 1:
              locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
              locale[iat][l][n][0] *= 0.5;
              locale[iat][l][n][1] += locale[iat][l][n][0];
					  	break;

					  case 2:
					  	for(int is=0; is<NSPIN; is++)
                locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
					  	break;

					  case 4: //SOC
					  	locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
					  	break;

					  default:
					  	cout << "Not supported NSPIN parameter" << endl;
					  	exit(0);			
					}

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
					locale_save[iat][l][n][0] = locale[iat][l][n][0];
					locale_save[iat][l][n][1] = locale[iat][l][n][1];	

          locale[iat][l][n][0].zero_out();
          locale[iat][l][n][1].zero_out();				
				}
			}			
		}
	}
	
	//=================Part 1======================
	//call PBLAS routine to calculate the product of the S and density matrix
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;

	vector<double> srho(ParaO.nloc);
	for(int is=0; is<NSPIN; is++)
	{
		// srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_gamma(iw,nu)
		pdgemm_(&transN, &transT,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				LM.Sloc, &one_int, &one_int, ParaO.desc, 
				LOC.wfc_dm_2d.dm_gamma.at(is).c, &one_int, &one_int, ParaO.desc,
				&beta,
				&srho[0], &one_int, &one_int, ParaO.desc);

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

	  								if( (nu>=0) && (mu>=0) )
	  								{																																																
	  									int m0_all = m0 + (2*l+1)*ipol0;
	  									int m1_all = m0 + (2*l+1)*ipol1;

	  									locale[iat][l][n][is](m0,m1) += srho[irc]/4.0;														
	  								}

	  								if( (nu_prime>=0) && (mu_prime>=0) )
	  								{
	  									int m0_all = m0 + (2*l+1)*ipol0;
	  									int m1_all = m0 + (2*l+1)*ipol1;
  
	  									locale[iat][l][n][is](m0,m1) += srho[irc_prime]/4.0;
	  								}
	  							}			
	  						}
	  					}
	  				}

	  				matrix temp(locale[iat][l][n][is]);
	  				MPI_Allreduce( &temp(0,0), &locale[iat][l][n][is](0,0), (2*l+1)*NPOL*(2*l+1)*NPOL,
	  								      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	  				// for the case spin independent calculation
	  				switch(NSPIN)
	  				{
	  				  case 1:	
	  				  	locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                locale[iat][l][n][0] *= 0.5;
                locale[iat][l][n][1] += locale[iat][l][n][0];
	  				  	break;

	  				  case 2:
	  				  	locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
	  				  	break;

	  				  default:
	  				  	cout << "Not supported NSPIN parameter" << endl;
	  				  	exit(0);			
	  				}

	  			}//end for(n)
	  			//this->print(it, iat, 2*L+1, N);
	  		}//L
	  	}//ia
	  }//it

	}//is

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
					
					if(NSPIN==1 || NSPIN==2)
					{
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
					}
					else if(NSPIN==4) //SOC
					{
						for(int m0=0; m0<2*l+1; m0++)
						{
							for(int ipol0=0; ipol0<NPOL; ipol0++)
							{
								const int m0_all = m0 + (2*l+1)*ipol0;

								for(int m1=0; m1<2*l+1; m1++)
								{
									for(int ipol1=0; ipol1<NPOL; ipol1++)
									{
										int m1_all = m1 + (2*l+1)*ipol1;
										ofdftu << fixed << setw(12) << setprecision(8) << locale.at(iat).at(l).at(n).at(0)(m0_all, m1_all);
									}
								}
								ofdftu << endl;
							}							
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
				cout << "DFTU::read_occup_m. Can not find the file initial_onsite.dm . Please check your initial_onsite.dm" << endl;
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
							
							if(NSPIN==1 || NSPIN==2)
							{
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
							else if(NSPIN==4) //SOC
							{
								double value = 0.0;
								for(int m0=0; m0<2*L+1; m0++)
								{
									for(int ipol0=0; ipol0<NPOL; ipol0++)
									{
										const int m0_all = m0 + (2*L+1)*ipol0;

										for(int m1=0; m1<2*L+1; m1++)
										{		
											for(int ipol1=0; ipol1<NPOL; ipol1++)	
											{		
												int m1_all = m1 + (2*L+1)*ipol1;	
												ifdftu >> value;
												locale.at(iat).at(L).at(zeta).at(0)(m0_all, m1_all) = value;		
											}				 	
										}
										ifdftu.ignore(150, '\n');
									}											
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
										
					if(NSPIN==1 || NSPIN==2)
					{
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
					else if(NSPIN==4) //SOC
					{
						for(int m0=0; m0<2*L+1; m0++)
						{
							for(int ipol0=0; ipol0<NPOL; ipol0++)
							{
								const int m0_all = m0 + (2*L+1)*ipol0;

								for(int m1=0; m1<2*L+1; m1++)
								{		
									for(int ipol1=0; ipol1<NPOL; ipol1++)	
									{		
										int m1_all = m1 + (2*L+1)*ipol1;	
										
										MPI_Bcast(&locale.at(iat).at(l).at(n).at(0)(m0_all, m1_all), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);		
									}				 	
								}
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
						if(NSPIN==1 || NSPIN==2)
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
						else if(NSPIN==4) //SOC
						{
							double nm_trace = 0.0;
							double nm2_trace = 0.0;

							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int ipol0=0; ipol0<NPOL; ipol0++)
								{
							 		const int m0_all = m0 + (2*l+1)*ipol0;
									nm_trace += this->locale.at(iat).at(l).at(n).at(0)(m0_all, m0_all);

									for(int m1=0; m1<2*l+1; m1++)
									{
										for(int ipol1=0; ipol1<NPOL; ipol1++)
										{
											int m1_all = m1 + (2*l+1)*ipol1;

											nm2_trace += this->locale.at(iat).at(l).at(n).at(0)(m0_all,m1_all)*this->locale.at(iat).at(l).at(n).at(0)(m1_all, m0_all);
										}
									}
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
						for(int ipol1=0; ipol1<NPOL; ipol1++)
						{
							const int m1_all = m1 + ipol1*(2*l+1);
							for(int m2=0; m2<2*l+1; m2++)
							{	
								for(int ipol2=0; ipol2<NPOL; ipol2++)
								{
									const int m2_all = m2 + ipol2*(2*l+1);

									if(NSPIN==1 || NSPIN==2)
									{
										for(int is=0; is<2; is++)
										{
											double VU = 0.0;
											VU = get_onebody_eff_pot(T, iat, l, n, is, m1_all, m2_all, cal_type, 0);
											EU_dc += VU*this->locale.at(iat).at(l).at(n).at(is)(m1_all,m2_all);
										}
									}	
									else if(NSPIN==4) //SOC
									{										
										double VU = 0.0;
										VU = get_onebody_eff_pot(T, iat, l, n, 0, m1_all, m2_all, cal_type, 0);
										EU_dc += VU*this->locale.at(iat).at(l).at(n).at(0)(m1_all,m2_all);										
									}
								}					
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

void DFTU::cal_eff_pot_mat_complex(const int ik, const int istep, complex<double>* eff_pot)
{
	TITLE("DFTU", "cal_eff_pot_mat");

 	if((CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax") && (!INPUT.omc) && istep==0 && this->iter_dftu==1) return;
	
	int spin = kv.isk[ik];

	ZEROS(eff_pot, ParaO.nloc);

	//ofs_running << "dftu.cpp "<< __LINE__  << endl;
	//=============================================================
	//   PART2: call pblas to calculate effective potential matrix
	//=============================================================
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
  const complex<double> alpha_c(1.0,0.0), beta_c(0.0,0.0), half_c(0.5,0.0), one_c(1.0,0.0);

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
			int m1_all = m1 + (2*L1+1)*ipol1;
			int m2_all = m2 + (2*L2+1)*ipol2;
			double val = get_onebody_eff_pot(T1, iat1, L1, n1, spin, m1_all, m2_all, cal_type, 1);
			VU.at(irc) = complex<double>(val, 0.0);
		}
	}

	pzgemm_(&transN, &transN,
		&NLOCAL, &NLOCAL, &NLOCAL,
		&half_c, 
		VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc,
		LM.Sloc2, &one_int, &one_int, ParaO.desc,
		&beta_c,
		eff_pot, &one_int, &one_int, ParaO.desc);

  for(int irc=0; irc<ParaO.nloc; irc++)
    VU[irc] = eff_pot[irc];
  
  // pztranc(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
  pztranc_(&NLOCAL, &NLOCAL, 
           &one_c, 
           &VU[0], &one_int, &one_int, ParaO.desc, 
           &one_c, 
           eff_pot, &one_int, &one_int, ParaO.desc);

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

void DFTU::cal_eff_pot_mat_real(const int ik, const int istep, double* eff_pot)
{
	TITLE("DFTU", "cal_eff_pot_mat");

 	if((CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax") && (!INPUT.omc) && istep==0 && this->iter_dftu==1) return;
	
	int spin = kv.isk[ik];

	ZEROS(eff_pot, ParaO.nloc);

	//ofs_running << "dftu.cpp "<< __LINE__  << endl;
	//=============================================================
	//   PART2: call pblas to calculate effective potential matrix
	//=============================================================
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0, half=0.5, one=1.0;

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
			VU.at(irc) = get_onebody_eff_pot(T1, iat1, L1, n1, spin, m1, m2, cal_type, 1);
		}
	}
		
	// The first term
	pdgemm_(&transN, &transN,
		&NLOCAL, &NLOCAL, &NLOCAL,
		&half, 
		VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc, 
		LM.Sloc, &one_int, &one_int, ParaO.desc,
		&beta,
		eff_pot, &one_int, &one_int, ParaO.desc);

  for(int irc=0; irc<ParaO.nloc; irc++)
    VU[irc] = eff_pot[irc];
  
  // pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
  pdtran_(&NLOCAL, &NLOCAL, 
          &one, 
          &VU[0], &one_int, &one_int, ParaO.desc, 
          &one, 
          eff_pot, &one_int, &one_int, ParaO.desc);

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
					
					if(NSPIN==1 || NSPIN==2)
					{
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
					}
					else if(NSPIN==4) //SOC
					{
						for(int m0=0; m0<2*l+1; m0++)
						{
							for(int ipol0=0; ipol0<NPOL; ipol0++)
							{
								const int m0_all = m0 + (2*l+1)*ipol0;

								for(int m1=0; m1<2*l+1; m1++)
								{
									for(int ipol1=0; ipol1<NPOL; ipol1++)
									{
										int m1_all = m1 + (2*l+1)*ipol1;
										ofs_running << fixed << setw(12) << setprecision(8) << locale.at(iat).at(l).at(n).at(0)(m0_all, m1_all);
									}
								}
								ofs_running << endl;
							}							
						}
					}					
				}//n
			}//l
		}//I
	}//T
	ofs_running << "//=======================================================//" << endl;
	
	return;
}

void DFTU::cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR)
{
  const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0, one=1.0, half=0.5;

  for(int i=0; i<ParaO.nloc; i++) HR[i] = 0.0;

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

			VU.at(irc) = get_onebody_eff_pot(T1, iat1, L1, n1, ispin, m1, m2, cal_type, 1);
	  }//ic
  }//ir

	// The first term
	pdgemm_(&transN, &transN,
		&NLOCAL, &NLOCAL, &NLOCAL,
		&half, 
		VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc, 
		SR, &one_int, &one_int, ParaO.desc,
		&beta,
		HR, &one_int, &one_int, ParaO.desc);

	for(int irc=0; irc<ParaO.nloc; irc++)
		VU[irc] = HR[irc];

  pdtran_(&NLOCAL, &NLOCAL, 
          &one, 
          &VU[0], &one_int, &one_int, ParaO.desc, 
          &one, 
          HR, &one_int, &one_int, ParaO.desc);

  return;
}

void DFTU::cal_eff_pot_mat_R_complex_double(
  const int ispin, complex<double>* SR, complex<double>* HR)
{
  const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const complex<double> alpha(1.0,0.0), beta(0.0,0.0);
  const complex<double> zero(0.0,0.0), half(0.5,0.0), one(1.0,0.0);

  for(int i=0; i<ParaO.nloc; i++) HR[i] = zero;

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

			int m1_all = m1 + (2*L1+1)*ipol1;
			int m2_all = m2 + (2*L2+1)*ipol2;

			double val = get_onebody_eff_pot(T1, iat1, L1, n1, ispin, m1_all, m2_all, cal_type, 1);

			VU.at(irc) = complex<double>(val, 0.0);
		}
	}
	// vector<complex<double>> potm_tmp(ParaO.nloc, complex<double>(0.0, 0.0));

	// The first term
	pzgemm_(&transN, &transN,
		&NLOCAL, &NLOCAL, &NLOCAL,
		&half, 
		VECTOR_TO_PTR(VU), &one_int, &one_int, ParaO.desc,
		SR, &one_int, &one_int, ParaO.desc,
		&beta,
		HR, &one_int, &one_int, ParaO.desc);

	for(int irc=0; irc<ParaO.nloc; irc++)
	  VU[irc] = HR[irc];

  pztranc_(&NLOCAL, &NLOCAL, 
           &one, 
           &VU[0], &one_int, &one_int, ParaO.desc, 
           &one, 
           HR, &one_int, &one_int, ParaO.desc);

  return;
}

void DFTU::folding_overlap_matrix(const int ik, complex<double>* Sk)
{
  TITLE("DFTU","folding_overlap_matrix"); 
	// timer::tick("DFTU","folding_overlap_matrix");

  ZEROS(Sk, ParaO.nloc);

	int iat = 0;
	int index = 0;
	Vector3<double> dtau;
	Vector3<double> tau1;
	Vector3<double> tau2;

	Vector3<double> dtau1;
	Vector3<double> dtau2;
	Vector3<double> tau0;

	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GridD.Find_atom(tau1);
			GridD.Find_atom(ucell, tau1, T1, I1);
			Atom* atom1 = &ucell.atoms[T1];
			const int start = ucell.itiaiw2iwt(T1,I1,0);

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];

				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

				bool adj = false;

				if(distance < rcut) 
				{
					adj = true;
				}
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0); 
						const int I0 = GridD.getNatom(ad0); 
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						double distance1 = dtau1.norm() * ucell.lat0;
						double distance2 = dtau2.norm() * ucell.lat0;

						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							adj = true;
							break;
						}
					}
				}

				if(adj) // mohan fix bug 2011-06-26, should not be '<='
				{
					// (3) calculate the nu of atom (T2, I2)
					const int start2 = ucell.itiaiw2iwt(T2,I2,0);
					//------------------------------------------------
					// exp(k dot dR)
					// dR is the index of box in Crystal coordinates
					//------------------------------------------------
					Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z); 
					const double arg = ( kv.kvec_d[ik] * dR ) * TWO_PI;
					//const double arg = ( kv.kvec_d[ik] * GridD.getBox(ad) ) * TWO_PI;
					const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );

					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					//--------------------------------------------------
					for(int ii=0; ii<atom1->nw*NPOL; ii++)
					{
						// the index of orbitals in this processor
						const int iw1_all = start + ii;
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<atom2->nw*NPOL; jj++)
						{
							int iw2_all = start2 + jj;
							const int nu = ParaO.trace_loc_col[iw2_all];

							if(nu<0)continue;
							//const int iic = mu*ParaO.ncol+nu;
              int iic;
              if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
              {
                  iic=mu+nu*ParaO.nrow;
              }
              else
              {
                  iic=mu*ParaO.ncol+nu;
              }

							//########################### EXPLAIN ###############################
							// 1. overlap matrix with k point
							// LM.SlocR = < phi_0i | phi_Rj >, where 0, R are the cell index
							// while i,j are the orbital index.

							// 2. H_fixed=T+Vnl matrix element with k point (if Vna is not used).
							// H_fixed=T+Vnl+Vna matrix element with k point (if Vna is used).
							// LM.Hloc_fixed = < phi_0i | H_fixed | phi_Rj>

							// 3. H(k) |psi(k)> = S(k) | psi(k)> 
							// Sloc2 is used to diagonalize for a give k point.
							// Hloc_fixed2 is used to diagonalize (eliminate index R).
							//###################################################################
							
							if(NSPIN!=4)
							{
								Sk[iic] += LM.SlocR[index] * kphase;
							}
							else
							{
								Sk[iic] += LM.SlocR_soc[index] * kphase;
							}
							++index;

						}//end jj
					}//end ii
				}
			}// end ad
			++iat;
		}// end I1
	} // end T1

	assert(index==LNNR.nnr);

  // timer::tick("DFTU","folding_overlap_matrix");
	return;
}