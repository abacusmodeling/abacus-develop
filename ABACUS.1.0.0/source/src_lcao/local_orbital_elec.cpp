#include "local_orbital_elec.h"
#include "diago_lcao_matrix.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
// mohan add on test 2010-01-13
//#include "../src_develop/src_onscaling/on_tests.h"
#include "update_input.h"
#include "src_pw/chi0_hilbert.h"
#include "lcao_vna.h"
#include "evolve_lcao_matrix.h"
#include "../src_pw/berryphase.h"
#include "../src_pw/toWannier90.h"

int Local_Orbital_Elec::iter = 0;
double Local_Orbital_Elec::avg_iter = 0.0; 
int Local_Orbital_Elec::istep = 0;

void Local_Orbital_Elec::scf(const int &istep)
{
	TITLE("Local_Orbtial_Elec","scf");
	timer::tick("Local_Orbital_Elec","scf",'D');

//	ofs_running << " Local_Orbital_Electrons_relaxation begin :"<<endl;
//	ofs_running << " Self-consistent Calculation "<<endl;

	// (1) calculate ewald energy.
	en.ewld = en.ewald();

    set_ethr(); //mohan add 2012-02-08

	Symmetry_rho srho;
	for(int is=0; is<NSPIN; is++)
	{
		srho.begin(is);
	}

//	cout << scientific;
//	cout << setiosflags(ios::fixed);

	if(OUT_LEVEL=="ie" ||OUT_LEVEL=="m") //xiaohui add "m" option, 2015-09-16
	{
		if(COLOUR && MY_RANK==0)
		{
			printf( " \e[33m%-7s\e[0m", "ITER");	
			printf( "\e[33m%-15s\e[0m", "ETOT(Ry)");	
			if(NSPIN==2)
			{
				printf( "\e[33m%-10s\e[0m", "TMAG");	
				printf( "\e[33m%-10s\e[0m", "AMAG");	
			}
			printf( "\e[33m%-14s\e[0m", "DRHO2");	
			printf( "\e[33m%-15s\e[0m", "ETOT(eV)");	
			//printf( "\e[33m%-11s\e[0m", "BAND(Ry)");	
			//printf( "\e[33m%-11s\e[0m", "H(Ry)");	
			//printf( "\e[33m%-11s\e[0m", "EXC(Ry)");	
			//if(DIAGO_TYPE=="cg")
			//{
			//	printf( "\e[33m%-11s\e[0m", "CG_ITER");	
			//}
			printf( "\e[33m%-11s\e[0m\n", "TIME(s)");	
		}
		else
		{	
			cout << " " << setw(7)<< "ITER";
			//cout << setw(15) << "ETOT(Ry)";

			if(NSPIN==2)
			{
				cout<<setw(10)<<"TMAG";
				cout<<setw(10)<<"AMAG";
			}

                        cout << setw(15) << "ETOT(eV)";
                        cout << setw(15) << "EDIFF(eV)";
                        cout << setw(11) << "DRHO2";
			//cout << setw(11) << "BAND(RY)";
			//cout << setw(11) << "H(RY)";
			//cout << setw(11) << "EXC(RY)"; 
			//if(DIAGO_TYPE=="cg")
			//{
			//	cout << setw(11) << "CG_ITER";
			//}
			//else if(DIAGO_TYPE=="selinv")
			//{
			//	cout << setw(11) << "SIAO_ITER";
			//}
			cout << setw(11) << "TIME(s)" << endl;
		}
	}

	for(Local_Orbital_Elec::iter=1; iter<=NITER; iter++)
	{
        if(CALCULATION=="scf")
        {
            ofs_running
            << "\n LCAO ALGORITHM ------------- ELEC=" << setw(4) << iter
            << "--------------------------------\n";

            ofs_warning
            << "\n LCAO ALGORITHM ------------- ELEC=" << setw(4) << iter
            << "--------------------------------\n";
        }
        else if(CALCULATION=="relax" || CALCULATION=="cell-relax")
		{
			ofs_running 
			<< "\n LCAO ALGORITHM ------------- ION=" << setw(4) << istep+1 
			<< "  ELEC=" << setw(4) << iter 
			<< "--------------------------------\n";

			ofs_warning 
			<< "\n LCAO ALGORITHM ------------- ION=" << setw(4) << istep+1 
			<< "  ELEC=" << setw(4) << iter 
			<< "--------------------------------\n";
		}
		else if(CALCULATION=="md")
		{
			ofs_running 
			<< "\n LCAO ALGORITHM ------------- MD=" << setw(4) << istep+1 
			<< "  ELEC=" << setw(4) << iter 
			<< "--------------------------------\n";

			ofs_warning 
			<< "\n LCAO ALGORITHM ------------- MD=" << setw(4) << istep+1 
			<< "  ELEC=" << setw(4) << iter 
			<< "--------------------------------\n";
		}

		time_t time_start, time_finish;
		clock_t clock_start;

		string ufile = "CHANGE";
		Update_input UI;
		UI.init(ufile);
			
		
		time_start= std::time(NULL);
		clock_start = std::clock();
		conv_elec = false;//mohan add 2008-05-25

		// mohan add 2010-07-16
		// used for pulay mixing.
		if(iter==1)	chr.set_new_e_iteration(true);
		else		chr.set_new_e_iteration(false);

		// set converged threshold, 
		// automatically updated during self consistency, only for CG.
        this->update_ethr(iter);
        if(FINAL_SCF && iter==1)
        {
            init_mixstep_final_scf();
            //chr.irstep=0;
            //chr.idstep=0;
            //chr.totstep=0;
        }

		// mohan update 2012-06-05
		en.calculate_harris(1);

		// mohan move it outside 2011-01-13
		// first need to calculate the weight according to 
		// electrons number.
		// mohan add iter > 1 on 2011-04-02
		// because the en.ekb has not value now.
		// so the smearing can not be done.
		if(iter>1)Occupy::calculate_weights();
		
		if(wf.start_wfc == "file")
		{
			if(iter==1)
			{
				cout << " WAVEFUN -> CHARGE " << endl;

				// The occupation should be read in together.
				// Occupy::calculate_weights(); //mohan add 2012-02-15
				
				// calculate the density matrix using read in wave functions
				// and the ncalculate the charge density on grid.
				LOC.sum_bands();
				// calculate the local potential(rho) again.
				// the grid integration will do in later grid integration.


				// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// a puzzle remains here.
				// if I don't renew potential,
				// The dr2 is very small.
				// OneElectron, Hartree and
				// Exc energy are all correct
				// except the band energy.
				//
				// solved by mohan 2010-09-10
				// there are there rho here:
				// rho1: formed by read in orbitals.
				// rho2: atomic rho, used to construct H
				// rho3: generated by after diagonalize
				// here converged because rho3 and rho1
				// are very close.
				// so be careful here, make sure
				// rho1 and rho2 are the same rho.
				// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				pot.v_of_rho(chr.rho, en.ehart, en.etxc, en.vtxc, pot.vr);
				en.delta_escf();
				if (vext == 0)	pot.set_vrs(pw.doublegrid);
				else		pot.set_vrs_tddft(pw.doublegrid, istep);
			}
		}

		//fuxiang add 2016-11-1
		if(tddft==1 && iter == 1)
		{
			this->WFC_init = new complex<double>**[kv.nks];
			for(int ik=0; ik<kv.nks; ik++)
			{
				this->WFC_init[ik] = new complex<double>*[NBANDS];
			}
			for(int ik=0; ik<kv.nks; ik++)
			{
				for(int ib=0; ib<NBANDS; ib++)
				{
					this->WFC_init[ik][ib] = new complex<double>[NLOCAL];
				}
			}
			if(istep>=1)
			{
				for (int ik=0; ik<kv.nks; ik++)
				{
					for (int ib=0; ib<NBANDS; ib++)
					{
						for (int i=0; i<NLOCAL; i++)
						{
							WFC_init[ik][ib][i] = LOWF.WFC_K[ik][ib][i];
						}
					}
				}
			}
			else
			{
				for (int ik=0; ik<kv.nks; ik++)
				{
					for (int ib=0; ib<NBANDS; ib++)
					{
						for (int i=0; i<NLOCAL; i++)
						{
							WFC_init[ik][ib][i] = complex<double>(0.0,0.0);
						}
					}
				}
			}
		}	
		
		// calculate exact-exchange
		switch(xcf.iexch_now)						// Peize Lin add 2018-10-30
		{
			case 5:    case 6:   case 9:
				if( !exx_global.info.separate_loop )				
				{
					exx_lcao.cal_exx_elec();
				}
				break;
		}		
		
		// (1) calculate the bands.
		cal_bands(istep);

//		for(int ib=0; ib<NBANDS; ++ib)
//		{
//			cout << ib+1 << " " << wf.ekb[0][ib] << endl;
//		}

		//-----------------------------------------------------------
		// only deal with charge density after both wavefunctions.
		// are calculated.
		//-----------------------------------------------------------
		if(GAMMA_ONLY_LOCAL && NSPIN == 2 && CURRENT_SPIN == 0)continue;


		if(conv_elec) 
		{
			timer::tick("Local_Orbital_Elec","scf",'D');
			return;
		}
		
		en.eband  = 0.0;
		en.ef     = 0.0;
		en.ef_up  = 0.0;
		en.ef_dw  = 0.0;

		// demet is included into eband.
		//if(DIAGO_TYPE!="selinv")
		{
			en.demet  = 0.0;
		}

		// (2)
		chr.save_rho_before_sum_band();
		
		// (3) sum bands to calculate charge density
		Occupy::calculate_weights();

		if (ocp == 1)
		{
			for (int ik=0; ik<kv.nks; ik++)
			{
				for (int ib=0; ib<NBANDS; ib++)
				{
					wf.wg(ik,ib)=ocp_kb[ik*NBANDS+ib];
				}
			}
		}

/*		if(istep==0)
		{
			for (int ib=0; ib < 3; ib++)
			{
				wf.wg(0,ib) = 2.0;
			}
			wf.wg(0,3) = 1.75;
			wf.wg(0,4) = 0.25;
		}
    		for (int ik = 0; ik < kv.nks; ik++)
    		{
        		for (int ib = 0; ib < NBANDS; ib++)
			{
				cout << "wf.wg(ik, ib): " << wf.wg(ik, ib) << endl;
			}
    		}
*/

		for(int ik=0; ik<kv.nks; ++ik)
		{
			en.print_band(ik);
		}

		// if selinv is used, we need this to calculate the charge
		// using density matrix.
		LOC.sum_bands();

		// add exx
		en.set_exx();		// Peize Lin add 2016-12-03

		// (4) mohan add 2010-06-24
		// using new charge density.
		en.calculate_harris(2);
		
		// (5) symmetrize the charge density
		Symmetry_rho srho;
		for(int is=0; is<NSPIN; is++)
		{
			srho.begin(is);
		}

		// (6) compute magnetization, only for spin==2
        mag.compute_magnetization();

		// resume codes!
		//-------------------------------------------------------------------------
		//this->LOWF.init_Cij( 0 ); // check the orthogonality of local orbital.
		//chr.sum_band(); use local orbital in plane wave basis to calculate bands.
		// but must has evc first!
		//-------------------------------------------------------------------------

		// (7) calculate delta energy
		en.deband = en.delta_e();

		// (8) Mix charge density
		chr.mix_rho(dr2,0,DRHO2,iter,conv_elec);

		// (9) Calculate new potential according to new Charge Density.
	
		if(conv_elec || iter==NITER)
		{ 
			if(pot.out_potential<0) //mohan add 2011-10-10
			{
				pot.out_potential = -2;
			}
		}

		if(!conv_elec)
		{
			// option 1
			pot.v_of_rho(chr.rho, en.ehart, en.etxc, en.vtxc, pot.vr);
			en.delta_escf();

			// option 2
			//------------------------------
			// mohan add 2012-06-08
			// use real E_tot functional.
			//------------------------------
			/*
			pot.v_of_rho(chr.rho_save, en.ehart, en.etxc, en.vtxc, pot.vr);
			en.calculate_etot();
			en.print_etot(conv_elec, istep, iter, dr2, 0.0, ETHR, avg_iter,0);
			pot.v_of_rho(chr.rho, en.ehart, en.etxc, en.vtxc, pot.vr);
			en.delta_escf();
			*/
		}
		else
		{
			pot.v_of_rho(chr.rho, en.ehart, en.etxc, en.vtxc, pot.vnew);
			//(used later for scf correction to the forces )
			pot.vnew -= pot.vr;
			en.descf = 0.0;
		}

		//-----------------------------------
		// output charge density for tmp
		//-----------------------------------
		for(int is=0; is<NSPIN; is++)
		{
			const int precision = 3;
			
			stringstream ssc;
			ssc << global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
			chr.write_rho( is, iter, ssc.str(), precision );//mohan add 2007-10-17

			stringstream ssd;

			if(GAMMA_ONLY_LOCAL)
			{
				ssd << global_out_dir << "tmp" << "_SPIN" << is + 1 << "_DM";
			}
			else
			{
				ssd << global_out_dir << "tmp" << "_SPIN" << is + 1 << "_DM_R";
			}
			LOC.write_dm( is, iter, ssd.str(), precision );

			stringstream ssp;
			ssp << global_out_dir << "tmp" << "_SPIN" << is + 1 << "_POT";
			pot.write_potential( is, iter, ssp.str(), pot.vr, precision );
		}

		// (10) add Vloc to Vhxc.
		if(vext == 0)	pot.set_vrs(pw.doublegrid);
		else		pot.set_vrs_tddft(pw.doublegrid, istep);
	
		time_finish=std::time(NULL);
		double duration = (double)(clock() - clock_start) / CLOCKS_PER_SEC;

		// (11) calculate the total energy.
		en.calculate_etot();
		en.print_etot(conv_elec, istep, iter, dr2, duration, ETHR, avg_iter);
	
		en.etot_old = en.etot;

		if (conv_elec || iter==NITER)
		{
			//--------------------------------------
			// output charge density for converged,
			// 0 means don't need to consider iter,
			//--------------------------------------
			if( chi0_hilbert.epsilon)                                    // pengfei 2016-11-23
			{
				cout <<"eta = "<<chi0_hilbert.eta<<endl;
				cout <<"domega = "<<chi0_hilbert.domega<<endl;
				cout <<"nomega = "<<chi0_hilbert.nomega<<endl;
				cout <<"dim = "<<chi0_hilbert.dim<<endl;
				//cout <<"oband = "<<chi0_hilbert.oband<<endl;
				chi0_hilbert.Chi();
			}

			for(int is=0; is<NSPIN; is++)
			{
				const int precision = 3;

        		stringstream ssc;
       			ssc << global_out_dir << "SPIN" << is + 1 << "_CHG";
        		chr.write_rho( is, 0, ssc.str() );//mohan add 2007-10-17

				stringstream ssd;
				if(GAMMA_ONLY_LOCAL)
				{
					ssd << global_out_dir << "SPIN" << is + 1 << "_DM";
				}
				else
				{
					ssd << global_out_dir << "SPIN" << is + 1 << "_DM_R";
				}
				LOC.write_dm( is, 0, ssd.str(), precision );

				stringstream ssp;
				ssp << global_out_dir << "SPIN" << is + 1 << "_POT";
				pot.write_potential( is, 0, ssp.str(), pot.vrs, precision );

				//fuxiang add 2017-03-15
				stringstream sse;
				sse << global_out_dir << "SPIN" << is + 1 << "_DIPOLE_ELEC";
				chr.write_rho_dipole( is, 0, sse.str());

			}
			
			iter_end(ofs_running);

			if(conv_elec)
			{
 				//xiaohui add "OUT_LEVEL", 2015-09-16
				if(OUT_LEVEL != "m") ofs_running << setprecision(16);
				if(OUT_LEVEL != "m") ofs_running << " EFERMI = " << en.ef * Ry_to_eV << " eV" << endl; 
				if(OUT_LEVEL=="ie")
				{
					ofs_running << " " << global_out_dir << " final etot is " << en.etot * Ry_to_eV << " eV" << endl; 
				}
			}
			else
			{
				ofs_running << " !! convergence has not been achieved @_@" << endl;
				if(OUT_LEVEL=="ie" || OUT_LEVEL=="m") //xiaohui add "m" option, 2015-09-16
				cout << " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!" << endl;
			}

//			DONE(ofs_running,"ELECTRONS CONVERGED!");
			timer::tick("Local_Orbital_Elec","scf",'D');
			return;
		}
	}

	if (tddft==1)
	{
		delete[] WFC_init;
	}

	timer::tick("Local_Orbital_Elec","scf",'D');
	return;
}

void Local_Orbital_Elec::nscf(void)
{
	TITLE("Local_Orbital_Elec","nscf");

	cout << " NON-SELF CONSISTENT CALCULATIONS" << endl;
	
	time_t time_start= std::time(NULL);
//	complex<double>*** WFC_init;
	this->WFC_init = new complex<double>**[kv.nks];
	for(int ik=0; ik<kv.nks; ik++)
	{
		this->WFC_init[ik] = new complex<double>*[NBANDS];
	}
	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int ib=0; ib<NBANDS; ib++)
		{
			this->WFC_init[ik][ib] = new complex<double>[NLOCAL];
		}
	}

	for (int ik=0; ik<kv.nks; ik++)
	{
		for (int ib=0; ib<NBANDS; ib++)
		{
			for (int i=0; i<NLOCAL; i++)
			{
				WFC_init[ik][ib][i] = complex<double>(0.0,0.0);
			}
		}
	}

	switch(exx_lcao.info.hybrid_type)			// Peize Lin add 2018-08-14
	{
		case Exx_Global::Hybrid_Type::HF:
		case Exx_Global::Hybrid_Type::PBE0:
		case Exx_Global::Hybrid_Type::HSE:
			exx_lcao.cal_exx_elec_nscf();
			break;
	}

//	cal_bands(istep, WFC_init);
	cal_bands(istep);
	time_t time_finish=std::time(NULL);
	OUT_TIME("cal_bands",time_start, time_finish);

    ofs_running << " end of band structure calculation " << endl;
    ofs_running << " band eigenvalue in this processor (eV) :" << endl;

    for (int ik = 0; ik < kv.nks; ik++)
    {
        if (NSPIN==2)
        {
            if (ik == 0) ofs_running << " spin up :" << endl;
            if (ik == ( kv.nks / 2)) ofs_running << " spin down :" << endl;
        }
        //out.printV3(ofs_running, kv.kvec_c[ik]);

		ofs_running << " k-points" << ik+1 << "(" << kv.nkstot << "): " << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z << endl;

        for (int ib = 0; ib < NBANDS; ib++)
        {			
            ofs_running << " spin" << kv.isk[ik]+1 << "final_state " << ib+1 << " " << wf.ekb[ik][ib] * Ry_to_eV << " " << wf.wg(ik, ib)*kv.nks << endl;
        //    cout << " spin" << kv.isk[ik]+1 << "final_state " << ib+1 << " " << wf.ekb[ik][ib] * Ry_to_eV << " " << wf.wg(ik, ib)*kv.nks << endl;
        }
		ofs_running << endl;
    }
	
	// add by jingan in 2018.11.7
	if(CALCULATION == "nscf" && INPUT.towannier90)
    {
        toWannier90 myWannier(kv.nkstot,ucell.G);
        myWannier.init_wannier();
    }
	
	// add by jingan
	if (BERRY_PHASE && SYMMETRY == 0)
    {
    	berryphase bp;
		bp.Macroscopic_polarization();
    }

	return;
}

#include "../src_parallel/subgrid_oper.h"
void Local_Orbital_Elec::cal_bands(const int &istep)
{
	TITLE("Local_Orbital_Elec","cal_bands");
	timer::tick("Local_Orbital_Elec","cal_bands",'E');

	if(GAMMA_ONLY_LOCAL)
	{
		assert(NSPIN == kv.nks);
	}
	else
	{
		int start_spin = -1;
		UHM.GK.reset_spin(start_spin);
		UHM.GK.allocate_pvpR();
	}

	ofs_running << " "  <<setw(8) << "K-point" << setw(15) << "Time(Sec)"<< endl;
	ofs_running << setprecision(6) << setiosflags(ios::fixed) << setiosflags(ios::showpoint);
	for(int ik=0; ik<kv.nks; ik++)
	{	
		//cout << " ik=" << ik+1 << " spin=" << kv.isk[ik] << endl;
		//out.printV3(kv.kvec_c[ik]);
		
		//-----------------------------------------
		//(1) prepare data for this k point.
		// copy the local potential from array.
		//-----------------------------------------
		if(NSPIN==2)CURRENT_SPIN = kv.isk[ik];
		wf.npw = kv.ngk[ik];
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			pot.vrs1[ir] = pot.vrs( CURRENT_SPIN, ir);
		}
		
		//--------------------------------------------
		//(2) check if we need to calculate 
		// pvpR = < phi0 | v(spin) | phiR> for
		// a new spin.
		//--------------------------------------------
		if(!GAMMA_ONLY_LOCAL)
		{
			if(CURRENT_SPIN == UHM.GK.get_spin() )
			{
				//ofs_running << " Same spin, same vlocal integration." << endl;
			}
			else
			{
				ofs_running << " (spin change)" << endl;
				UHM.GK.reset_spin( CURRENT_SPIN );
			
				// if you change the place of the following code,
				// rememeber to delete the #include	
				if(VL_IN_H)
				{
					if(VNA==0)
					{
						// vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
						UHM.GK.cal_vlocal_k(pot.vrs1,GridT);
						if(NONCOLIN) //&& DOMAG)
						{//added by zhengdy-soc, for non-collinear case
							for(int is=1;is<4;is++)
							{
								for(int ir=0; ir<pw.nrxx; ir++)
								{
									pot.vrs1[ir] = pot.vrs( is, ir);
								}
								UHM.GK.cal_vlocal_k(pot.vrs1, GridT, is);
							}
						}
					}
					else if(VNA>1)
					{
						LCAO_Vna lv;
						lv.smooth_vl2();
					}
					else if(VNA==1)
					{
						LCAO_Vna::smooth_vl1();
					}
				}
			}
		}


		// some preparation
		clock_t start=clock();
//		cout << " init s matrix = " << UHM.init_s << endl;
		if(!UHM.init_s)
    	{
    	    WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		clock_t start_nscf=clock();

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		if(GAMMA_ONLY_LOCAL)
		{
			UHM.calculate_Hgamma(ik);						// Peize Lin add ik 2016-12-03

			// SGO: sub_grid_operation
			SGO.cal_totwfc();


			//--------------------------------------
			// DIAG GROUP OPERATION HERE
			//--------------------------------------
			if(DCOLOR==0)
			{
				/*
				//xiaohui modify 2013-09-02
				//if(LINEAR_SCALING==1)
				//{
				//	Diago_LCAO_Matrix DLM;
				//	// write the wave functions into LOWF.WFC_GAMMA.


				//	// mohan add 2011-04-15
				//	if(BFIELD)
				//	{
				//		//cout << " solve comlex matrix()" << endl;
				//		DLM.solve_complex_matrix(ik, SGO.totwfc_B[0], istep);
				//	}
				//	else
				//	{
				//		//DLM.solve_double_matrix(ik, LOWF.WFC_GAMMA[CURRENT_SPIN]);
				//		// the temperary array totwfc only have one spin direction.
				//		DLM.solve_double_matrix(ik, SGO.totwfc[0]);
				//	}
				//}
				//else if(LINEAR_SCALING==2)
				//{	
				//	if(NURSE)
				//	{
				//		// the full matrix method.
				//		//ON.run_gamma(ik);
				//	}
				//	else
				//	{
				//		//ON.run_gamma_sparse(ik);
				//	}
				//}
				//else
				//{
				//	WARNING_QUIT("Local_Orbital_Elec::cal_bands","check LINEAR_SCALING");
				//} xiaohui modify 2013-09-02. Attention...
				*/
				
				//xiaohui add 2013-09-02
				Diago_LCAO_Matrix DLM;

				if(BFIELD)
				{
					//cout << " solve comlex matrix()" << endl;
					DLM.solve_complex_matrix(ik, SGO.totwfc_B[0], LOC.wfc_dm_2d.wfc_k[ik]);
				}
				else
				{
					//DLM.solve_double_matrix(ik, LOWF.WFC_GAMMA[CURRENT_SPIN]);
					// the temperary array totwfc only have one spin direction.
					DLM.solve_double_matrix(ik, SGO.totwfc[0], LOC.wfc_dm_2d.wfc_gamma[ik]);
				} //xiaohui add 2013-09-02. Attention...
			}
			else
			{
#ifdef __MPI
				ofs_running << " no diagonalization." << endl;
#else
				cout << " DCOLOR=" << DCOLOR << endl;
				WARNING_QUIT("Local_Orbital_Elec::cal_bands","no diagonalization");
#endif

			}
#ifdef __MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			// distribute the wave functions again.
			if(!BFIELD)
			{
				SGO.dis_subwfc();
			}
			else
			{
				SGO.dis_subwfc_complex();
			}
		}//end gamma
		// with k points
		else
		{
			timer::tick("Efficience","each_k");
			timer::tick("Efficience","H_k");
			UHM.calculate_Hk(ik);
			timer::tick("Efficience","H_k");

			// write the wave functions into LOWF.WFC_K[ik].
			bool diago = true;
			if (tddft == 1 && istep >= 1) diago = false;
			if(diago)
			{
				timer::tick("Efficience","diago_k");
				Diago_LCAO_Matrix DLM;
				DLM.solve_complex_matrix(ik, LOWF.WFC_K[ik], LOC.wfc_dm_2d.wfc_k[ik]);
				timer::tick("Efficience","diago_k");
			}
			else
			{
				timer::tick("Efficience","evolve_k");
				Evolve_LCAO_Matrix ELM;
				ELM.evolve_complex_matrix(ik, LOWF.WFC_K[ik], WFC_init[ik]);
				timer::tick("Efficience","evolve_k");
			}
			timer::tick("Efficience","each_k");
		}


//		cout << " In Local_Orbital_Elec" << endl;
//		for(int ie=0; ie<NBANDS; ie++)
//		{
//			cout << " e[" << ie << "]=" << wf.ekb[ik][ie] * Ry_to_eV<< endl;
//		}

		//BLOCK_HERE("Bands after diagonalization");

		/* 
		 * test : use Cij to and evc to calculate new evc
		 bool use_blas = true;
		 if(use_blas)
		 {
		// do evc(NBANDS, npwx) = ONE * Cij(NBANDS, NLOCAL) * wf.wanf2[0](NLOCAL, npwx) + ZERO
		zgemm('N', // op( Cij ) = Cij
		'N',  // op( wf.wanf2[0] ) = wf.wanf2[0]
		NBANDS, // row of evc
		kv.ngk[0], // column of evc
		NLOCAL, // between Cij and wf.wanf2[0]
		ONE, // see formula
		this->LOWF.c,
		NBANDS, // row of op( Cij ), if 'T' or 'C', in fact this means the column number
		wf.wanf2[0],
		NLOCAL, // row of op( wf.wanf2[0] )
		ZERO, // see formula
		wf.evc[0],
		NBANDS // first dimension of evc
		);
		}
		else
		{
		wf.evc[0] = this->LOWF.c * wf.wanf2[0];
		}
		 */

		clock_t finish_nscf=clock();
		const double duration = static_cast<double>(finish_nscf - start_nscf) / CLOCKS_PER_SEC;


		ofs_running << setw(8) << ik+1 << setw(15) << duration << endl;

//		ofs_running << " TIME FOR K=" << ik << ": " << duration << endl;
	}
			
	// delete at last time.
	if(!GAMMA_ONLY_LOCAL)
	{
		UHM.GK.destroy_pvpR();
	}
	timer::tick("Local_Orbital_Elec","cal_bands",'E');
	return;	
}

void Local_Orbital_Elec::init_mixstep_final_scf(void)
{
    TITLE("Local_Orbital_Elec","init_mixstep_final_scf");

    chr.irstep=0;
    chr.idstep=0;
    chr.totstep=0;

    return;
}
