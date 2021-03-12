#include "LOOP_ions.h"
#include "src_pw/global.h"
#include "src_parallel/parallel_orbitals.h"
#include "src_pdiag/pdiag_double.h"
#include "LCAO_nnr.h"
#include "FORCE_STRESS.h"
#include "src_global/global_function.h"
#include "src_io/hs_matrix.h"
#include "src_io/cal_r_overlap_R.h"
#include "src_ions/variable_cell.h" // mohan add 2021-02-01
#include "src_ri/exx_abfs.h"
#include "src_ri/exx_opt_orb.h"
#include "ELEC_scf.h"
#include "src_global/sltk_atom_arrange.h"
#include "src_pw/vdwd2.h"

LOOP_ions::LOOP_ions()
{}

LOOP_ions::~LOOP_ions() 
{}

void LOOP_ions::opt_ions(void)
{
    TITLE("LOOP_ions","opt_ions"); 
    timer::tick("LOOP_ions","opt_ions",'B'); 
		
    if(OUT_LEVEL=="i")
    {
        cout << setprecision(12);
        cout<< " " << setw(7)<< "ISTEP"
        <<setw(5)<< "NE"
        <<setw(18)<< "ETOT(eV)"
        <<setw(10)<< "dE(meV)"
        <<setw(10)<< "F(eV/A)"
        <<setw(10)<< "T(MIN)"
        <<endl;
    }

    // Geometry optimization algorithm setup.
    if(FORCE)
    {
        //Ions_Move_Methods 
        IMM.allocate();
        //Charge_Extrapolation
        CE.allocate_ions();
    }

	// pengfei Li 2018-05-14
    if(STRESS)
    {
		// allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }

    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;
    while(istep <= NSTEP && !stop)
    {
        time_t estart = time(NULL);
	
		// xiaohui add "m" option, 2015-09-16
        if(OUT_LEVEL=="ie" || OUT_LEVEL=="m")
        {
            cout << " ---------------------------------------------------------" << endl;
            if(CALCULATION=="relax") 
            {
                cout << " RELAX IONS : " << istep << endl;
            }
            else if(CALCULATION=="cell-relax")
            {
                cout << " RELAX CELL : " << stress_step << endl;
                cout << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
            }
            else if(CALCULATION=="scf")
            {
                cout << " SELF-CONSISTENT : " << endl;
            }
            cout << " ---------------------------------------------------------" << endl;

            ofs_running << " ---------------------------------------------------------" << endl;
            if(CALCULATION=="relax")
            {
                ofs_running << " RELAX IONS : " << istep << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
            else if(CALCULATION=="cell-relax")
            {
                ofs_running << " RELAX CELL : " << stress_step << endl;
                ofs_running << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
            else if(CALCULATION=="scf")
            {
                ofs_running << " SELF-CONSISTENT" << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
        }

		// solve electronic structures in terms of LCAO
		// mohan add 2021-02-09
		LOE.solve_elec_stru(this->istep);

		
		time_t eend = time(NULL);

        //xiaohui add 2014-07-07, for second-order extrapolation
        int iat=0;
        if(CALCULATION=="relax" || CALCULATION=="cell-relax")
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

        if(pot.out_potential == 2)
        {
            stringstream ssp;
            stringstream ssp_ave;
            ssp << global_out_dir << "ElecStaticPot";
            ssp_ave << global_out_dir << "ElecStaticPot_AVE";
            pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        if(ParaO.out_hsR) 
		{
			this->output_HS_R(); //LiuXh add 2019-07-15
		}

        time_t fstart = time(NULL);
        if (CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            stop = this->force_stress(istep, force_step, stress_step);
        }            
        time_t fend = time(NULL);


        //xiaohui add 2014-07-07, for second-order extrapolation
        iat=0;
        if(FORCE)
        {
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
        }
		
        if(OUT_LEVEL=="i")
        {
            double etime_min = difftime(eend, estart)/60.0;
            double ftime_min = difftime(fend, fstart)/60.0;
            stringstream ss;
            ss << MOVE_IONS << istep;

            cout << setiosflags(ios::scientific)
            << " " << setw(7) << ss.str()
            << setw(5) << ELEC_scf::iter
            << setw(18) << setprecision(6) << en.etot * Ry_to_eV;

            cout << setprecision(2) << setiosflags(ios::scientific)
            << setw(10) << IMM.get_ediff() * Ry_to_eV * 1000
            << setw(10) << IMM.get_largest_grad() * Ry_to_eV / BOHR_TO_A;
            //<< setw(12) << IMM.get_trust_radius();

            cout << resetiosflags(ios::scientific)
//            << setw(8) << IMM.get_update_iter()
            << setprecision(2) << setw(10) << etime_min + ftime_min;
            cout << endl;
        }

//#ifdef __MPI //2015-09-06, xiaohui
	//2015-05-07, 2015-10-01
        //atom_arrange::delete_vector( SEARCH_RADIUS );
//#endif //2015-09-06, xiaohui

//2015-09-16
//#ifdef __MPI
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i=0;i<ucell.ntype;i++)
//    {
//        ucell.atoms[i].bcast_atom(); // bcast tau array
//    }
//#endif

        ++istep;
    }

    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
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
        		Force_LCAO FL;
                matrix stress_lcao;
                stress_lcao.create(3,3);
                FL.cal_stress(stress_lcao);
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

	// mohan update 2021-02-10
    hm.orb_con.clear_after_ions();

    timer::tick("LOOP_ions","opt_ions",'B'); 
    return;
}

//bool LOOP_ions::force_stress(void)
bool LOOP_ions::force_stress(const int &istep, int &force_step, int &stress_step)
{
    TITLE("LOOP_ions","force_stress");
    if(!FORCE && !STRESS)
    {
        return 1;
    }
    timer::tick("LOOP_ions","force_stress",'D');
	matrix fcs;
	matrix scs;
	Force_Stress_LCAO FSL;
	FSL.allocate (); 
	FSL.getForceStress(FORCE, STRESS, TEST_FORCE, TEST_STRESS, fcs, scs);
	//--------------------------------------------------
	// only forces are needed, no stresses are needed
	//--------------------------------------------------
    if(FORCE && !STRESS)
    {

#ifdef __MPI //2015-10-01, xiaohui
        atom_arrange::delete_vector( SEARCH_RADIUS );
#endif //2015-10-01, xiaohui

        if(CALCULATION=="relax") 
        {
            IMM.cal_movement(istep, istep, fcs, en.etot);

            if(IMM.get_converged() || (istep==NSTEP))
            {
                return 1; // 1 means converged
            }
            else // ions are not converged
            {
                CE.istep = istep;
                CE.extrapolate_charge();

                if(pot.extra_pot=="dm")
                {
                }
                else
                {
                    pot.init_pot( istep );
                }
            }
            return 0;
        }
        else
        {
            return 1;
        }

        // mohan update 2013-04-11
        // setup the structure factor
        // and do the density extraploation.
        // for both ionic iteratoin and
        // force calculations.

        //xiaohui modify 2014-08-09
        //pw.setup_structure_factor();

        // charge extrapolation if istep>0.
        //xiaohui modify 2014-08-09
        //CE.extrapolate_charge();

/*xiaohui modify 2014-08-09
        if(pot.extra_pot==4)
        {
            // done after grid technique.
        }
        else
        {
            pot.init_pot( istep );
        }
xiaohui modify 2014-08-09*/
    }

    static bool converged_force = false;
    static bool converged_stress = false;

    if(!FORCE&&STRESS)
    {

#ifdef __MPI
		atom_arrange::delete_vector( SEARCH_RADIUS );
#endif
		if(CALCULATION=="cell-relax")
		{
           	LCM.cal_lattice_change(stress_step, scs, en.etot);
           	converged_stress = LCM.get_converged();
           	if(converged_stress)
           	{
               	return 1;
           	}
           	else
           	{
               	Variable_Cell::init_after_vc();
               	pot.init_pot(stress_step);

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

//#ifdef __MPI
        atom_arrange::delete_vector( SEARCH_RADIUS );
//#endif
        
        //if(CALCULATION=="relax") IMM.cal_movement(istep, FL.fcs, en.etot);
        if(CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            IMM.cal_movement(istep, force_step, fcs, en.etot);

            if(IMM.get_converged())
            {
                force_step = 1;


			    if(CALCULATION=="cell-relax")
			    {
            	    LCM.cal_lattice_change(stress_step, scs, en.etot);
            	    converged_stress = LCM.get_converged();
            	    if(converged_stress)
            	    {
                	    return 1;
            	    }
            	    else
            	    {
                	    Variable_Cell::init_after_vc();
                	    pot.init_pot(stress_step);

                	    ++stress_step;
                	    return 0;
                    }
                }
                else
                {
                    return 1;
                }

#ifdef __MPI
            //atom_arrange::delete_vector( SEARCH_RADIUS );
#endif
            }
            else
            {
#ifdef __MPI
            //atom_arrange::delete_vector( SEARCH_RADIUS );
#endif
                //CE.istep = istep;
                CE.istep = force_step;
                CE.extrapolate_charge();

                if(pot.extra_pot=="dm")//xiaohui modify 2015-02-01
                {
                    // done after grid technique.
                }
                else
                {
                    pot.init_pot( istep );
                }
                ++force_step;
                return 0;
            }
        }
        else
        {
            return 1;
        }
    }

    return 0;

    timer::tick("LOOP_ions","force_stress",'D');
}

void LOOP_ions::final_scf(void)
{
    TITLE("LOOP_ions","final_scf"); 

    FINAL_SCF = true;
    Variable_Cell::final_calculation_after_vc();
    atom_arrange::set_sr_NL();
    atom_arrange::search( SEARCH_RADIUS );
    GridT.set_pbc_grid(
        pw.ncx, pw.ncy, pw.ncz,
        pw.bx, pw.by, pw.bz,
        pw.nbx, pw.nby, pw.nbz,
        pw.nbxx, pw.nbzp_start, pw.nbzp);

    // (2) If k point is used here, allocate HlocR after atom_arrange.
    if(!GAMMA_ONLY_LOCAL)
    {
        // For each atom, calculate the adjacent atoms in different cells 
        // and allocate the space for H(R) and S(R).
        LNNR.cal_nnr();
        LM.allocate_HS_R(LNNR.nnr);
        
		// need to first calculae lgd.
        // using GridT.init.
        LNNR.cal_nnrg(GridT);
    }

    // (4) set the augmented orbitals index.
    // after ParaO and GridT, 
    // this information is used to calculate
    // the force.
    LOWF.set_trace_aug(GridT);
		
    // (5) init density kernel
    // (6) init wave functions.
    if(GAMMA_ONLY_LOCAL)
    {
        // here we reset the density matrix dimension.
        LOC.allocate_gamma(GridT);
    }
    else
    {
        LOWF.allocate_k(GridT);
        LOC.allocate_DM_k();
    }

    UHM.set_lcao_matrices();

    if(vdwd2_para.flag_vdwd2)							//Peize Lin add 2014-04-04, update 2021-03-09
    {
        Vdwd2 vdwd2(ucell,vdwd2_para);
        vdwd2.cal_energy();
        en.evdw = vdwd2.get_energy();
    }
	else if(vdwd3.vdwD3)							//jiyy add 2019-05-18
    {
        vdwd3.energy();
        en.evdw = vdwd3.energy_result;
    }											  
    
	ELEC_scf es;
	es.scf(0);

    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
    {
        ofs_running << "\n\n --------------------------------------------" << endl;
        ofs_running << setprecision(16);
        ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
        ofs_running << " --------------------------------------------\n\n" << endl;
    }

    return;
}

void LOOP_ions::output_HS_R(void)
{
    TITLE("LOOP_ions","output_HS_R"); 
    timer::tick("LOOP_ions","output_HS_R",'D'); 
	
	// add by jingan for out r_R matrix 2019.8.14
	if(INPUT.out_r_matrix)
	{
		cal_r_overlap_R r_matrix;
		r_matrix.init();
		r_matrix.out_r_overlap_R(NSPIN);
	}

    if(NSPIN==1||NSPIN==4)
    {
        UHM.calculate_STN_R();
        UHM.GK.cal_vlocal_R(0);
        UHM.GK.distribute_pvpR_tr();
        HS_Matrix::save_HSR_tr(0);
    }
    ///*
    else if(NSPIN==2)
    {
        UHM.calculate_STN_R();
        for(int ik=0; ik<kv.nks; ik++)
        {
            if(ik==0 || ik==kv.nks/2)
            {
                if(NSPIN==2)CURRENT_SPIN = kv.isk[ik];
                for(int ir=0; ir<pw.nrxx; ir++)
                {
                    pot.vrs1[ir] = pot.vrs( CURRENT_SPIN, ir);
                }
        	    	
                if(!GAMMA_ONLY_LOCAL)
                {
                    if(VL_IN_H)
                    {
						//UHM.GK.cal_vlocal_k(pot.vrs1,GridT);
						UHM.GK.cal_vlocal_k(pot.vrs1,GridT,CURRENT_SPIN);
                    }
                }
                UHM.GK.cal_vlocal_R(CURRENT_SPIN);
                UHM.GK.distribute_pvpR_tr();
                HS_Matrix::save_HSR_tr(CURRENT_SPIN);
            }
        }
    }
    //*/
    if(!GAMMA_ONLY_LOCAL) //LiuXh 20181011
    {
        UHM.GK.destroy_pvpR();
    } //LiuXh 20181011

    timer::tick("LOOP_ions","output_HS_R",'D'); 
    return;
}
