#include "tools.h"
#include "ions.h"
#include "forces.h"
#include "stress.h"
#include "algorithms.h"
#include "../src_pw/global.h" // use chr.
#include "vdwd2.h"
#include "vdwd3.h"

#include "../src_pw/pw_complement.h"
#include "../src_pw/pw_basis.h"

void Ions::opt_ions_pw(void)
{
	TITLE("Ions","opt_ions_pw");
	timer::tick("Ions","opt_ions_pw",'C');
	
	if(OUT_LEVEL=="i")
	{
		cout << setprecision(12);
    	cout<< " " << setw(7)<< "ISTEP" 
		<<setw(5)<< "NE"
		<<setw(15)<< "ETOT(eV)"
		<<setw(15)<< "EDIFF(eV)"
        <<setw(15)<< "MAX_F(eV/A)"
        <<setw(15)<< "TRADIUS(Bohr)"
		<<setw(8)<< "UPDATE"
		<<setw(11)<< "ETIME(MIN)"
		<<setw(11)<< "FTIME(MIN)"
        <<endl;
	}

	// allocation for ion movement.	
	if(FORCE)
	{
		IMM.allocate();
		CE.allocate();
	}

	if(STRESS)                    // pengfei Li 2018-05-14
	{
		LCM.allocate();
	}

    this->istep = 1;
	int force_step = 1;           // pengfei Li 2018-05-14
	int stress_step = 1;
	bool stop= false;
	
    while(istep <= NSTEP && !stop)
    {
		time_t estart = time(NULL);

		if(OUT_LEVEL=="ie")
		{
	
        	cout << " -------------------------------------------" << endl;
			if(CALCULATION=="relax") //pengfei 2014-10-13
			{
        		cout << " STEP OF ION RELAXATION : " << istep << endl;
			}
            else if(CALCULATION=="cell-relax")
            {
                cout << " RELAX CELL : " << stress_step << endl;
                cout << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
                cout << " ---------------------------------------------------------" << endl;
            }
			else if(CALCULATION=="scf") //add 4 lines 2015-09-06, xiaohui
			{
        			cout << " SELF-CONSISTENT : " << endl;
			}
			else if(CALCULATION=="md")
			{
        		cout << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
			}
        	cout << " -------------------------------------------" << endl;

        	ofs_running << " -------------------------------------------" << endl;
			if(CALCULATION=="relax")
			{
        		ofs_running << " STEP OF ION RELAXATION : " << istep << endl;
			}
            else if(CALCULATION=="cell-relax")
            {
                ofs_running << " RELAX CELL : " << stress_step << endl;
                ofs_running << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
			else if(CALCULATION=="md")
			{
        		ofs_running << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
			}
        	ofs_running << " -------------------------------------------" << endl;
		}

		
		
		
		if(vdwd2.vdwD2)											//Peize Lin add 2014-04-03, update 2019-04-26
		{
			vdwd2.energy();
		}
		if(vdwd3.vdwD3)											//jiyy add 2019-05-18
		{
			vdwd3.energy();
		}																										 
		
        if (CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax" || CALCULATION=="cell-relax")  // pengfei 2014-10-13
        {
			if( Exx_Global::Hybrid_Type::No==exx_global.info.hybrid_type  )
			{			
				this->self_consistent(istep-1);
			}
			else if( Exx_Global::Hybrid_Type::Generate_Matrix == exx_global.info.hybrid_type )
			{
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			}
			else	// Peize Lin add 2019-03-09
			{
				if( exx_global.info.separate_loop )
				{
					for( size_t hybrid_step=0; hybrid_step!=exx_global.info.hybrid_step; ++hybrid_step )
					{
						this->self_consistent(istep-1);
						if( electrons::iter==1 || hybrid_step==exx_global.info.hybrid_step-1 )		// exx converge
							break;
						exx_global.info.set_xcfunc(xcf);							
						exx_lip.cal_exx();
					}						
				}
				else
				{
					this->self_consistent(istep-1);	
					exx_global.info.set_xcfunc(xcf);
					this->self_consistent(istep-1);
				}
			}
        }
        else if(CALCULATION=="nscf")
        {
            this->non_self_consistent();
        }
	

    int iat=0; //LiuXh add 20180619
    if(CALCULATION=="relax"|| CALCULATION=="md" || CALCULATION=="cell-relax")
    {
        for(int it = 0;it < ucell.ntype;it++)
        {
            Atom* atom = &ucell.atoms[it];
    	for(int ia =0;ia< ucell.atoms[it].na;ia++)
    	{
    	    CE.pos_old2[3*iat  ] = CE.pos_old1[3*iat  ];
    	    CE.pos_old2[3*iat+1] = CE.pos_old1[3*iat+1];
    	    CE.pos_old2[3*iat+2] = CE.pos_old1[3*iat+2];
    	    
    	    CE.pos_old1[3*iat  ] = CE.pos_now[3*iat  ];
    	    CE.pos_old1[3*iat+1] = CE.pos_now[3*iat+1];
    	    CE.pos_old1[3*iat+2] = CE.pos_now[3*iat+2];
    
    	    CE.pos_now[3*iat  ] = atom->tau[ia].x*ucell.lat0;
    	    CE.pos_now[3*iat+1] = atom->tau[ia].y*ucell.lat0;
    	    CE.pos_now[3*iat+2] = atom->tau[ia].z*ucell.lat0;
    
    	    iat++;
            }
        }
    }
	

		
		time_t eend = time(NULL);
		time_t fstart = time(NULL);
		//stop = this->force_stress(istep);
        if (CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
			stop = this->force_stress(istep, force_step, stress_step);    // pengfei Li 2018-05-14
		}
		time_t fend = time(NULL);

		if(OUT_LEVEL=="i")
		{
			double etime_min = difftime(eend, estart)/60.0; 
			double ftime_min = difftime(fend, fstart)/60.0; 
			stringstream ss;
			ss << MOVE_IONS << istep;
			
			cout << " " << setw(7) << ss.str() 
			<< setw(5) << this->iter 
			<< setw(15) << setprecision(6) << en.etot * Ry_to_eV 
			<< setw(15) << IMM.get_ediff() * Ry_to_eV
			<< setprecision(3)
			<< setw(15) << IMM.get_largest_grad() * Ry_to_eV / 0.529177
			<< setw(15) << IMM.get_trust_radius()
			<< setw(8) << IMM.get_update_iter()
			<< setprecision(2) << setw(11) << etime_min
			<< setw(11) << ftime_min << endl;
		}

		++istep;

    }

    if(CALCULATION=="scf" || CALCULATION=="relax")
    {
        ofs_running << "\n\n --------------------------------------------" << endl;
        ofs_running << setprecision(16);
        ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
        ofs_running << " --------------------------------------------\n\n" << endl;

/*
        if(STRESS)
        {
            if(stress_step==1)
            {
		Stress ss;
                ss.cal_stress();
                matrix stress;
                stress.create(3,3);

		PRESSURE = (ss.sigmatot[0][0]+ss.sigmatot[1][1]+ss.sigmatot[2][2])/3;
            }
            double pressure = PRESSURE;
            en.etot = en.etot + ucell.omega * pressure;

            ofs_running << "\n\n --------------------------------------------" << endl;
            ofs_running << setprecision(16);
            ofs_running << " !FINAL_ETOT_IS (+ P*V) " << en.etot * Ry_to_eV << " eV" << endl; 
            ofs_running << " --------------------------------------------\n\n" << endl;
        }
*/
    }

/*
    if(stop && STRESS) //LiuXh add 20180619
    {
        FINAL_SCF = true;
        Run_Frag::final_calculation_after_vc();
        this->self_consistent(0);
    }
*/

	if(OUT_LEVEL=="i")
	{
		cout << " ION DYNAMICS FINISHED :)" << endl;
	}

	timer::tick("Ions","opt_ions_pw",'C');
    return;
}

//bool Ions::force_stress(const int &istep)
bool Ions::force_stress(const int &istep, int &force_step, int &stress_step)  // pengfei Li 2018-05-14
{
	TITLE("Ions","force_stress");
	if(!FORCE && !STRESS)
	{
		return 1;
	}

	//if(STRESS)
	//{
		//calculate the stress
		//Stress ss;
		//ss.cal_stress();
		//change the latvec
	//}	

	if(FORCE&&!STRESS)
	{
		// (1) calculate the force.
		Forces fcs(ucell.nat);
		fcs.init();

		// (2) move the ions.
		bool converged = false;
		//if(CALCULATION=="md") //mohan add
		//{
		//	md.init_md(istep, fcs.force);
		//}
		//else	
		//{
			//IMM.cal_movement(istep, fcs.force, en.etot);
		//}
		if(CALCULATION=="relax")
		{
			IMM.cal_movement(istep, istep, fcs.force, en.etot);
			converged = IMM.get_converged();

            if(converged || (istep==NSTEP) ) 
            {
                return 1;
            }
            else
            {
                ofs_running << " Setup the structure factor in plane wave basis." << endl;
                pw.setup_structure_factor();

                ofs_running << " Setup the extrapolated charge." << endl;
                // charge extrapolation if istep>0.
                CE.extrapolate_charge();
			
                ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << endl;
                // calculate the new potential accordint to
                // the new charge density.
                pot.init_pot( istep );

                ofs_running << " Setup the new wave functions?" << endl;
                // newd() not needed now(if Q in r space, needed).
                wf.wfcinit();
                // mp_bcast
            }
        }
        else
        {
            return 1;
        }
    }

	static bool converged_force = false;              // pengfe Li  2018-05-14
	static bool converged_stress = false;

	if(!FORCE&&STRESS)
	{
		Stress ss;
		ss.cal_stress();
		matrix stress;
		stress.create(3,3);

		double unit_transform = 0.0;
		unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * eps8;
		double external_stress[3] = {PRESS1,PRESS2,PRESS3};
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				stress(i,j) = ss.sigmatot[i][j];
				//OUT(ofs_running,"stress(i,j)", stress(i,j)); //LiuXh modify 20180619
			}
			stress(i,i) = ss.sigmatot[i][i] - external_stress[i]/unit_transform;
		}
		PRESSURE = (ss.sigmatot[0][0]+ss.sigmatot[1][1]+ss.sigmatot[2][2])/3;
		if(CALCULATION=="cell-relax")
		{
			LCM.cal_lattice_change(stress_step, stress, en.etot);
			converged_stress = LCM.get_converged();
            //cout <<"converged_stress = "<<converged_stress<<endl;
            if(converged_stress)
            {
                return 1;
            }
            else
            {
                Run_Frag::frag_init_after_vc();
                //pot.init_pot(0);
                pot.init_pot(stress_step); //LiuXh add 20180619
                ofs_running << " Setup the new wave functions?" << endl; //LiuXh add 20180619
                wf.wfcinit(); //LiuXh add 20180619
                ++stress_step;
                return 0;
            }
        }
        else
        {
            return 1;
        }
    }

	if(FORCE&&STRESS)
	{
		//cout<<" istep  force_step  stress_step  converged_force  converged_stress = "<<istep<<"  "<<force_step<<"  "<<stress_step<<"  "<<converged_force<<"  "<<converged_stress<<endl;

		Forces fcs(ucell.nat);
		fcs.init();
		//IMM.cal_movement(force_step, fcs.force, en.etot);
        if(CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            IMM.cal_movement(istep, force_step, fcs.force, en.etot);
            converged_force = IMM.get_converged();

            //cout<<"converged_force = "<<converged_force<<endl;
            if(converged_force)
            {
                force_step = 1;

                Stress ss;
                ss.cal_stress();
                matrix stress;
                stress.create(3,3);

                double unit_transform = 0.0;
                unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * eps8;
                double external_stress[3] = {PRESS1,PRESS2,PRESS3};
                for(int i=0;i<3;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        stress(i,j) = ss.sigmatot[i][j];
                        //OUT(ofs_running,"stress(i,j)", stress(i,j)); //LiuXh modify 20180619
                    }
                    stress(i,i) = ss.sigmatot[i][i] - external_stress[i]/unit_transform;
                }
                PRESSURE = (ss.sigmatot[0][0]+ss.sigmatot[1][1]+ss.sigmatot[2][2])/3;

                if(CALCULATION=="cell-relax")
                {
                    LCM.cal_lattice_change(stress_step, stress, en.etot);
                    converged_stress = LCM.get_converged();
                    //cout <<"converged_stress = "<<converged_stress<<endl;
                    if(converged_stress)
                    {
                        return 1;
                    }
                    else
                    {
                        Run_Frag::frag_init_after_vc();
                        //pot.init_pot(0);
                        pot.init_pot(stress_step); //LiuXh add 20180619

                        ofs_running << " Setup the new wave functions?" << endl; //LiuXh add 20180619
                        wf.wfcinit(); //LiuXh add 20180619

                        ++stress_step;
                        return 0;
                    }
                }
                else
                {
                    return 1;
                }
            }
            else
            {
                //stress_step = 1;
                pw.setup_structure_factor();
                int iat=0; //LiuXh add 20180619
                for(int it = 0;it < ucell.ntype;it++)
                {
                    Atom* atom = &ucell.atoms[it];
                    for(int ia =0;ia< ucell.atoms[it].na;ia++)
                    {
                        CE.pos_next[3*iat  ] = atom->tau[ia].x*ucell.lat0;
                        CE.pos_next[3*iat+1] = atom->tau[ia].y*ucell.lat0;
                        CE.pos_next[3*iat+2] = atom->tau[ia].z*ucell.lat0;

                        iat++;
                    }
                }
                CE.istep = force_step;

                CE.extrapolate_charge();
                pot.init_pot( istep );
                wf.wfcinit();
                ++force_step;
                return 0;
            }
        }
        else
        {
            Stress ss;
            ss.cal_stress();
            matrix stress;
            stress.create(3,3);

            double unit_transform = 0.0;
            unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * eps8;
            double external_stress[3] = {PRESS1,PRESS2,PRESS3};
            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    stress(i,j) = ss.sigmatot[i][j];
                    //OUT(ofs_running,"stress(i,j)", stress(i,j)); //LiuXh modify 20180619
                }
                stress(i,i) = ss.sigmatot[i][i] - external_stress[i]/unit_transform;
            }
            PRESSURE = (ss.sigmatot[0][0]+ss.sigmatot[1][1]+ss.sigmatot[2][2])/3;
            return 1;
        }
    }

    return 0;
}


void Ions::extrapolate_wfcs()
{
    TITLE("Ions","extrapolate_wfcs");
    // wave function extrapolation:
    // wfc_order = 0 nothing is done
    // wfc_order = 2 first order extrapolation:
    // |psi(t+dt)> = 2*|psi(t)> - |psi(t-dt)>
    // wfc_order = 3 second order extrapolation:
    // |psi(t+dt)> = |psi(t)> +
    // + alpha0*( |psi(t)> - |psi(t-dt)> )
    // + beta0* ( |psi(t-dt)> - |psi(t-2*dt)> )


    return;
}
