#include "run_md_lcao.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_orbitals.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_lcao/LCAO_nnr.h"
#include "../src_lcao/FORCE_STRESS.h"
#include "../module_base/global_function.h"
#include "../src_io/write_HS.h"
#include "../src_io/cal_r_overlap_R.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../src_lcao/ELEC_scf.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"

Run_MD_LCAO::Run_MD_LCAO()
{}

Run_MD_LCAO::~Run_MD_LCAO() 
{}


void Run_MD_LCAO::opt_cell(void)
{
	TITLE("Run_MD_LCAO","opt_cell");

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    GlobalC::wf.allocate_ekb_wg(GlobalC::kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP
    GlobalC::UFFT.allocate();

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step=0;
    GlobalC::pot.init_pot(ion_step, GlobalC::pw.strucFac);

	
	opt_ions();
	return;
}


void Run_MD_LCAO::opt_ions(void)
{
    TITLE("Run_MD_LCAO","opt_ions"); 
    timer::tick("Run_MD_LCAO","opt_ions"); 
		
    if(GlobalV::OUT_LEVEL=="i")
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
    if(GlobalV::FORCE)
    {
        //Ions_Move_Methods 
        IMM.allocate();
        //Charge_Extrapolation
        CE.allocate_ions();
    }

	// pengfei Li 2018-05-14
    if(GlobalV::STRESS)
    {
		// allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }

    MD_basic mdb(INPUT.mdp, GlobalC::ucell);
    int mdtype = INPUT.mdp.mdtype;

    this->istep = 1;
    bool stop = false;
    while(istep <= GlobalV::NSTEP && !stop)
    {
        time_t estart = time(NULL);
	
		// xiaohui add "m" option, 2015-09-16
        if(GlobalV::OUT_LEVEL=="ie" || GlobalV::OUT_LEVEL=="m")
        {
			cout << " ---------------------------------------------------------" << endl;
			cout<<" Molecular Dynamics STEP "<< mdb.getRealStep()<<endl;
			cout << " ---------------------------------------------------------" << endl;
        }
		//----------------------------------------------------------
		// about vdw, jiyy add vdwd3 and linpz add vdwd2
		//----------------------------------------------------------	
		if(INPUT.vdw_method=="d2")
		{
			// setup vdwd2 parameters
			GlobalC::vdwd2_para.initial_parameters(INPUT);
	        GlobalC::vdwd2_para.initset(GlobalC::ucell);
        }
        if(INPUT.vdw_method=="d3_0" || INPUT.vdw_method=="d3_bj")
        {
            GlobalC::vdwd3_para.initial_parameters(INPUT);
        }
        // Peize Lin add 2014.04.04, update 2021.03.09
        if(GlobalC::vdwd2_para.flag_vdwd2)
        {
            Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
            vdwd2.cal_energy();
            GlobalC::en.evdw = vdwd2.get_energy();
        }
        // jiyy add 2019-05-18, update 2021.05.02
        else if(GlobalC::vdwd3_para.flag_vdwd3)
        {
            Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
            vdwd3.cal_energy();
            GlobalC::en.evdw = vdwd3.get_energy();
        }



		// solve electronic structures in terms of LCAO
		// mohan add 2021-02-09
		LOE.solve_elec_stru(this->istep);

		time_t eend = time(NULL);

        //xiaohui add 2014-07-07, for second-order extrapolation
		CE.update_all_pos(GlobalC::ucell);

		if(mdtype==1||mdtype==2)   
		{
			mdb.runNVT(istep);
		}
		else if(mdtype==0)  
		{
			mdb.runNVE(istep);
		}
        else if(mdtype==-1)
        {
            stop = mdb.runFIRE(istep);
        }
        else
        {
            WARNING_QUIT("opt_ions", "mdtype should be -1~2!");
        }

        if(GlobalC::pot.out_potential == 2)
        {
            stringstream ssp;
            stringstream ssp_ave;
            ssp << GlobalV::global_out_dir << "ElecStaticPot";
            ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
            GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        time_t fstart = time(NULL);
        time_t fend = time(NULL);

        //xiaohui add 2014-07-07, for second-order extrapolation
		CE.save_pos_next(GlobalC::ucell);

		//xiaohui add CE.istep = istep 2014-07-07
		CE.update_istep(istep);

		// charge extrapolation if istep>0.
		CE.extrapolate_charge();

		if(GlobalC::pot.extra_pot=="dm")//xiaohui modify 2015-02-01
		{
			// done after grid technique.
		}
		else
		{
			GlobalC::pot.init_pot( istep, GlobalC::pw.strucFac );
		}


        if(GlobalV::OUT_LEVEL=="i")
        {
            double etime_min = difftime(eend, estart)/60.0;
            double ftime_min = difftime(fend, fstart)/60.0;
            stringstream ss;
            ss << GlobalV::MOVE_IONS << istep;

            cout << setiosflags(ios::scientific)
            << " " << setw(7) << ss.str()
            << setw(5) << ELEC_scf::iter
            << setw(18) << setprecision(6) << GlobalC::en.etot * Ry_to_eV;

            cout << setprecision(2) << setiosflags(ios::scientific)
            << setw(10) << IMM.get_ediff() * Ry_to_eV * 1000
            << setw(10) << IMM.get_largest_grad() * Ry_to_eV / BOHR_TO_A;

            cout << resetiosflags(ios::scientific)
            << setprecision(2) << setw(10) << etime_min + ftime_min;
            cout << endl;
        }

        ++istep;
    }

    if(istep>1) 
	{
		final_scf();
	}

	// mohan update 2021-02-10
    GlobalC::LOWF.orb_con.clear_after_ions(UOT, ORB, INPUT.out_descriptor);

    timer::tick("Run_MD_LCAO","opt_ions"); 
    return;
}

void Run_MD_LCAO::final_scf(void)
{
    TITLE("Run_MD_LCAO","final_scf"); 

    GlobalV::FINAL_SCF = true;

    Variable_Cell::final_calculation_after_vc();

    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
		GlobalV::ofs_running, 
		GlobalV::OUT_LEVEL, 
		ORB.get_rcutmax_Phi(), 
		ORB.get_rcutmax_Beta(), 
		GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(
		GlobalV::SEARCH_PBC,
		GlobalV::ofs_running,
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS, 
		GlobalV::test_atom_input);

    GridT.set_pbc_grid(
        GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz,
        GlobalC::pw.bx, GlobalC::pw.by, GlobalC::pw.bz,
        GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbz,
        GlobalC::pw.nbxx, GlobalC::pw.nbzp_start, GlobalC::pw.nbzp);

    // (2) If k point is used here, allocate HlocR after atom_arrange.
    if(!GlobalV::GAMMA_ONLY_LOCAL)
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
    GlobalC::LOWF.set_trace_aug(GridT);
		
    // (5) init density kernel
    // (6) init wave functions.
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        // here we reset the density matrix dimension.
        GlobalC::LOC.allocate_gamma(GridT);
    }
    else
    {
        GlobalC::LOWF.allocate_k(GridT);
        GlobalC::LOC.allocate_DM_k();
    }

    UHM.set_lcao_matrices();

    if(GlobalC::vdwd2_para.flag_vdwd2) //Peize Lin add 2014-04-04, update 2021-03-09
    {
        Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
        vdwd2.cal_energy();
        GlobalC::en.evdw = vdwd2.get_energy();
    }
	else if(GlobalC::vdwd3_para.flag_vdwd3) //jiyy add 2019-05-18, update 2021-05-02
    {
        Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
        vdwd3.cal_energy();
        GlobalC::en.evdw = vdwd3.get_energy();
    }												  
    
	ELEC_scf es;
	es.scf(0);

    GlobalV::ofs_running << "\n\n --------------------------------------------" << endl;
    GlobalV::ofs_running << setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * Ry_to_eV << " eV" << endl; 
    GlobalV::ofs_running << " --------------------------------------------\n\n" << endl;

    return;
}
