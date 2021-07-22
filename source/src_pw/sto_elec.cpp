#include "tools.h"
#include "global.h"
#include "sto_elec.h" 
#include "occupy.h"
#include "../module_symmetry/symmetry_rho.h"
#include "../src_io/wf_io.h"
#include "H_Ewald_pw.h"

double Stochastic_Elec::avg_iter = 0;

Stochastic_Elec::Stochastic_Elec()
{
}

Stochastic_Elec::~Stochastic_Elec()
{
}


void Stochastic_Elec::scf_stochastic(const int &istep)
{
	timer::tick("Elec_Stochastic","scf_stochastic");

	// mohan update 2021-02-25
	H_Ewald_pw::compute_ewald(ucell,pw); 

    set_ethr();
    
	if(OUT_LEVEL=="ie")
	{
		cout << setprecision(12);
    	cout<< " " << setw(7)<< "ITER"; 

		if(NSPIN==2)
		{
			cout<<setw(10)<<"TMAG";
			cout<<setw(10)<<"AMAG";
		}
	
        cout<<setw(15)<< "ETOT(eV)"<<setw(15)<< "EDIFF(eV)"<<setw(11)<< "DRHO2"; // pengfei Li added 2015-1-31
		if(KS_SOLVER=="cg") 
		{
			cout<<setw(11)<<"CG_ITER";
		}

		cout<<setw(11)<< "TIME(S)";

        cout<<endl;
	}

	Symmetry_rho srho;
	for(int is=0; is<NSPIN; is++)
	{
		srho.begin(is, CHR, pw, Pgrid, symm);
	}

	// conv_elec is a member of Threshold_Elec
	this->conv_elec = false;

	clock_t start,finish;
	double duration = 0.0;	
	STO_WF.init();
	int chetype = 1;

	for (this->iter = 1;iter <= NITER;iter++)
    {
		ofs_running 
		<< "\n PW-STOCHASTIC ALGO --------- ION=" << setw(4) << istep + 1
		<< "  ELEC=" << setw(4) << iter 
		<< "--------------------------------\n";
		if(iter==1) 
		{
			CHR.new_e_iteration = true;
        }
		else 
		{
			CHR.new_e_iteration = false;
		}
		if(FINAL_SCF && iter==1)
        {
            CHR.irstep=0;
            CHR.idstep=0;
            CHR.totstep=0;
        }
		
		// record the start time.
        start=std::clock();
		
		
		//(1) set converged threshold, 
		// automatically updated during self consistency.
        //this->update_ethr(iter);
        if(FINAL_SCF && iter==1) 
		{
			ETHR = 1.0e-4/CHR.nelec; //smaller ETHR than KS-DFT
		}
        else 
		{
			if (iter == 2)
        	{
            	ETHR = 1.0e-4/CHR.nelec;
        	}
			ETHR = std::min( ETHR, 0.1*dr2/ std::max(1.0, CHR.nelec));
        }

        

		// mohan move harris functional to here, 2012-06-05
		// use 'rho(in)' and 'v_h and v_xc'(in)
		//en.calculate_harris(1);
	
		// first_iter_again:					// Peize Lin delete 2019-05-01
		
		//(2) calculate band energy using cg or davidson method.
		// output the new eigenvalues and wave functions.
		if(NBANDS > 0 && MY_POOL == 0)
		{
			this->c_bands(istep+1);
		}
		else
		{
			hm.hpw.init_k(0); //only GAMMA
			//In fact, hm.hpw.init_k has been done in wf.wfcinit();
		}
		kv.wk[0] = 2;// GAMMA temporary
		
        if (check_stop_now()) return;
        
		en.eband  = 0.0;
        en.demet  = 0.0;
        en.ef     = 0.0;
        en.ef_up  = 0.0;
        en.ef_dw  = 0.0;

		//(3) save change density as previous charge,
		// prepared fox mixing.
		if(MY_POOL == 0)
		{
        	CHR.save_rho_before_sum_band();
		}

		//prepare wavefunction&eband for other pools
#ifdef __MPI
		if(NBANDS > 0)
		{
			MPI_Bcast(wf.evc[0].c, wf.npwx*NBANDS*2, MPI_DOUBLE , 0, PARAPW_WORLD);
			MPI_Bcast(wf.ekb[0], NBANDS, MPI_DOUBLE, 0, PARAPW_WORLD);
		}
#endif
		
		//(4) calculate fermi energy.
		int ndim;
    	if(STO_WF.stotype == "pw")
		{
    	    ndim = wf.npw;
		}
    	else
		{
    	    ndim = pw.nrxx;
		}

		if(iter == 1)	
		{
			stoiter.init( ndim, chetype );
		}

		stoiter.orthog();
		stoiter.checkemm(iter);	//check and reset emax & emin
		stoiter.test();  //only for test
		stoiter.itermu(iter);


		//(5) calculate new charge density 
		// calculate KS rho.
		if(NBANDS > 0)
		{
			if(MY_POOL == 0)
			{
				CHR.sum_band();
			}
			else
			{
				for(int is=0; is<NSPIN; is++)
				{
					ZEROS(CHR.rho[is], pw.nrxx);
				}
			}
			//for(int is = 0; is < NSPIN; ++is)
			//{
			//	MPI_Bcast(CHR.rho[is], pw.nrxx, MPI_DOUBLE , 0,PARAPW_WORLD);
			//}
#ifdef __MPI
			MPI_Bcast(&en.eband,1, MPI_DOUBLE, 0,PARAPW_WORLD);
#endif
		}
		else
		{
			for(int is=0; is<NSPIN; is++)
			{
				ZEROS(CHR.rho[is], pw.nrxx);
			}
		}

	// calculate stochastic rho
		stoiter.sum_stoband();
		

		//(6) calculate the delta_harris energy 
		// according to new charge density.
		// mohan add 2009-01-23
		//en.calculate_harris(2);


		Symmetry_rho srho;
		for(int is=0; is<NSPIN; is++)
		{
			srho.begin(is, CHR, pw, Pgrid, symm);
		}


		//(8) deband is calculated from "output" charge density calculated 
		// in sum_band
		// need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
        en.deband = en.delta_e();
	
	
		// tr2_min used only in first scf iteraton
		double diago_error = 0.0;
		if(iter==1) 
		{
			// if 'dr2 < ETHR * nelec' happen,
			// in other word, 'dr2 < diago_error'
			// we update ETHR.
			diago_error = ETHR*std::max(1.0, CHR.nelec);
		}

		// if converged is achieved, or the self-consistent error(dr2)
		// is samller than the estimated error due to diagonalization(diago_error)
		// rhoin and rho are unchanged:
		// rhoin contain the input charge density and 
		// rho contain the output charge density.
		// in other cases rhoin contains the mixed charge density
		// (the new input density) while rho is unchanged.
		if(MY_POOL == 0)
		{
			CHR.mix_rho(dr2,diago_error,DRHO2,iter,conv_elec);
		}

#ifdef __MPI
		MPI_Bcast(&dr2, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		MPI_Bcast(&conv_elec, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		MPI_Bcast(CHR.rho[0], pw.nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

		//			if(MY_RANK==0)
		//			{
		//				ofs_mix << setw(5) << iter << setw(20) << dr2 << endl; 
		//			}

		if (iter==1)
		{
			if (dr2 < diago_error)
			{
				ofs_running << " Notice: Threshold on eigenvalues was too large.\n";

				WARNING("scf","Threshold on eigenvalues was too large.");
				ofs_running << " dr2=" << dr2 << " < diago_error=" << diago_error << endl;

				// update ETHR.
				ofs_running << " Origin ETHR = " << ETHR << endl;
				ETHR = dr2 / CHR.nelec;
				ofs_running << " New    ETHR = " << ETHR << endl;
				//                  goto first_iter_again;
			}
		}

        if (!conv_elec)
        {
			// not converged yet, calculate new potential from mixed charge density
            pot.vr = pot.v_of_rho(CHR.rho, CHR.rho_core);

			// because <T+V(ionic)> = <eband+deband> are calculated after sum
			// band, using output charge density.
			// but E_Hartree(en.ehart) and Exc(en.etxc) are calculated in v_of_rho above,
			// using the mixed charge density.
			// so delta_escf corrects for this difference at first order. 
            en.delta_escf();
        }
        else
        {
			for(int is=0; is<NSPIN; ++is)
			{
				for(int ir=0; ir<pw.nrxx; ++ir)
				{
					pot.vnew(is,ir) = pot.vr(is,ir);
				}
			}

			// the new potential V(PL)+V(H)+V(xc)
            pot.vr = pot.v_of_rho(CHR.rho, CHR.rho_core);

            //( vnew used later for scf correction to the forces )
	    	pot.vnew = pot.vr - pot.vnew;
            en.descf = 0.0;

        }
		
            

		// output for tmp.
		for(int is=0; is<NSPIN; is++)
		{
			stringstream ssc;
			ssc << global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
			CHR.write_rho(CHR.rho_save[is], is, iter, ssc.str(), 3);//mohan add 2007-10-17
		}
        
		

			pot.set_vr_eff();
        
        //print_eigenvalue(ofs_running);
		en.calculate_etot();

	
        finish=clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;

		en.print_etot(conv_elec, istep, iter, dr2, duration, ETHR, avg_iter);	
        if (conv_elec || iter==NITER)
        {
			for(int is=0; is<NSPIN; is++)
			{
        		stringstream ssc;
        		ssc << global_out_dir << "SPIN" << is + 1 << "_CHG";
        		CHR.write_rho(CHR.rho_save[is], is, 0, ssc.str() );//mohan add 2007-10-17
			}
              	 
			if(conv_elec)
			{
				//ofs_running << " convergence is achieved" << endl;			
				//ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
				ofs_running << " charge density convergence is achieved" << endl;
                                ofs_running << " final etot is " << en.etot * Ry_to_eV << " eV" << endl;
			}
			else
			{
				ofs_running << " convergence has NOT been achieved!" << endl;			
			}
							
			iter_end(ofs_running);
			timer::tick("Elec_Stochastic","scf_stochastic");
            return;
        }
		
        //if ( imix >= 0 )  CHR.rho = CHR.rho_save;
        //ofs_running << "\n start next iterate for idum ";
		
    } 
	
	timer::tick("Elec_Stochastic","scf_stochastic");
    return;
} // end electrons


bool Stochastic_Elec::check_stop_now(void)
{
    bool check_stop_now = false;

    if (check_stop_now)
    {
        conv_elec = false;
    }

    return check_stop_now;
} 


void Stochastic_Elec::c_bands(const int &istep)
{
	if (test_elec) TITLE("electrons","c_bands");
	timer::tick("Elec_Stochastic","c_bands");

	int precondition_type = 2;

	double *h_diag = new double[wf.npwx * NPOL];//added by zhengdy-soc
	ZEROS(h_diag, wf.npwx * NPOL);
    
	avg_iter = 0.0;
       
	ofs_running << " "  <<setw(8) << "K-point" << setw(15) << "CG iter num" << setw(15) << "Time(Sec)"<< endl;
	ofs_running << setprecision(6) << setiosflags(ios::fixed) << setiosflags(ios::showpoint);

	for (int ik = 0;ik < kv.nks;ik++)
	{
		hm.hpw.init_k(ik);

        //===========================================
        // Conjugate-Gradient diagonalization
        // h_diag is the precondition matrix
        // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
        //===========================================
        if (precondition_type==1)
        {
            for (int ig = 0;ig < wf.npw; ig++)
			{
				h_diag[ig] = max(1.0, wf.g2kin[ig]);
				if(NPOL==2) h_diag[ig+wf.npwx] = h_diag[ig];
			}
        }
        else if (precondition_type==2)
        {
            for (int ig = 0;ig < wf.npw; ig++)
			{
				h_diag[ig] = 1 + wf.g2kin[ig] + sqrt( 1 + (wf.g2kin[ig] - 1) * (wf.g2kin[ig] - 1));
				if(NPOL==2) h_diag[ig+wf.npwx] = h_diag[ig];
			}
        }
		//h_diag can't be zero!  //zhengdy-soc
		if(NPOL==2)
		{
			for(int ig = wf.npw;ig < wf.npwx; ig++)
			{
				h_diag[ig] = 1.0;
				h_diag[ig+ wf.npwx] = 1.0;
			}
		}

		clock_t start=clock();

		double avg_iter_k = 0.0;
		hm.diagH_pw(istep, this->iter, ik, h_diag, avg_iter_k);

		avg_iter += avg_iter_k;

		en.print_band(ik);

		clock_t finish=clock();
		const double duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;


		ofs_running << " " << setw(8) 
			<< ik+1 << setw(15) 
			<< avg_iter_k << setw(15) << duration << endl;
	}//End K Loop


	if(BASIS_TYPE=="pw")
	{
		//		ofs_running << " avg_iteri " << avg_iter << endl;
		//Parallel_Reduce::reduce_double_allpool(avg_iter); //mohan fix bug 2012-06-05
		//		ofs_running << " avg_iter_after " << avg_iter << endl;
		//avg_iter /= static_cast<double>(kv.nkstot);
	}

	delete [] h_diag;

	timer::tick("Elec_Stochastic","c_bands");

	return;
} 

