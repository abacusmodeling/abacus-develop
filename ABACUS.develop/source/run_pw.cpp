#include "run_pw.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "src_io/numerical_basis.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/print_info.h"
#include "src_pw/symmetry.h"

Run_pw::Run_pw(){}
Run_pw::~Run_pw(){}

void Run_pw::plane_wave_line(void)
{
    TITLE("Run_pw","plane_wave_line");
	timer::tick("Run_pw","plane_wave_line",'B');

    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put NLOCAL and NBANDS as input parameters
    ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running);
    //ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running, NLOCAL, NBANDS);
    DONE(ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (Symmetry::symm_flag)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    kv.set( symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info PI;
    PI.setup_parameters();

    // Initalize the plane wave basis set
    pw.gen_pw(ofs_running, ucell, kv);
    DONE(ofs_running,"INIT PLANEWAVE");
    cout << " UNIFORM GRID DIM     : " << pw.nx <<" * " << pw.ny <<" * "<< pw.nz << endl;
    cout << " UNIFORM GRID DIM(BIG): " << pw.nbx <<" * " << pw.nby <<" * "<< pw.nbz << endl;

    // mohan add 2010-10-10, just to test the symmetry of a variety
    // of systems.
    if(CALCULATION == "test")
    {
        Cal_Test::test_memory();
        QUIT();
    }

    // mohan add 2010-09-13
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    Pgrid.init(pw.ncx, pw.ncy, pw.ncz, pw.nczp,
        pw.nrxx, pw.nbz, pw.bz); // mohan add 2010-07-22, update 2011-05-04

//----------------------------------------------------------
// 1 read in initial data:
//   a lattice structure:atom_species,atom_positions,lattice vector
//   b k_points
//   c pseudopotential
// 2 setup planeware basis, FFT,structure factor, ...
// 3 initialize local and nonlocal pseudopotential in G_space
// 4 initialize charge desity and warefunctios in G_space
//----------------------------------------------------------

    //=====================================
    // init charge/potential/wave functions
    //=====================================
    CHR.allocate(NSPIN, pw.nrxx, pw.ngmc);
    pot.allocate(pw.nrxx);
  
    wf.allocate(kv.nks);

   
	UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    ppcell.init(ucell.ntype);

    //=====================
    // init hamiltonian
    //=====================
    hm.hpw.init(wf.npwx, NPOL, ppcell.nkb, pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc();
    DONE(ofs_running,"LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    ppcell.init_vnl();
    DONE(ofs_running,"NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    pot.init_pot(0);//atomic_rho, v_of_rho, set_vrs

    pot.newd();//once
    DONE(ofs_running,"INIT POTENTIAL");

    //==================================================
    // create ppcell.tab_at , for trial wave functions.
    //==================================================
    wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if ( NBANDS != 0 || (CALCULATION!="scf-sto" && CALCULATION!="relax-sto" && CALCULATION!="md-sto") )//qianrui add 
	{
    	wf.wfcinit();
    }
    
    

	switch(exx_global.info.hybrid_type)				// Peize Lin add 2019-03-09
	{
		case Exx_Global::Hybrid_Type::HF:
		case Exx_Global::Hybrid_Type::PBE0:
		case Exx_Global::Hybrid_Type::HSE:
			exx_lip.init(&kv, &wf, &pw, &UFFT, &ucell);
			break;
		case Exx_Global::Hybrid_Type::No:
			break;
		case Exx_Global::Hybrid_Type::Generate_Matrix:
		default:
			throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
	}

    DONE(ofs_running,"INIT BASIS");

	// ion optimization begins
	// electron density optimization is included in ion optimization
    
    Ions ions;
    ions.opt_ions_pw();



    // caoyu add 2020-11-24, mohan updat 2021-01-03
    if(BASIS_TYPE=="pw" && INPUT.out_descriptor==1)
    {
        Numerical_Descriptor nc;
        nc.output_descriptor(wf.evc, INPUT.lmax_descriptor);
        DONE(ofs_running,"GENERATE DESCRIPTOR FOR DEEPKS");
    }


    if(BASIS_TYPE=="pw" && winput::out_spillage) //xiaohui add 2013-09-01
    {
        //cout << "\n Output Spillage Information : " << endl;
        // calculate spillage value.
        if ( winput::out_spillage == 3)
        {
            // control here!!
            //LOCAL_BASIS = 0; xiaohui modify 2013-09-01
            BASIS_TYPE="pw"; //xiaohui add 2013-09-01
            //LOCAL_BASIS = 0;
            cout << " NLOCAL = " << NLOCAL << endl;

            for (int ik=0; ik<kv.nks; ik++)
            {
                wf.wanf2[ik].create(NLOCAL, wf.npwx);
                //if ( LOCAL_BASIS == 3 || LOCAL_BASIS == 0 ) xiaohui modify 2013-09-01
				if(BASIS_TYPE=="pw") //xiaohui add 2013-09-01. Attention! "LOCAL_BASIS==3"???
                {
					cout << " ik=" << ik + 1 << endl;

                    // mohan modify 2010-1-10
                    //LOCAL_BASIS=4; xiaohui modify 2013-09-1
                    BASIS_TYPE="lcao_in_pw"; //xiaohui add 2013-09-01. Attention! How about "BASIS_TYPE=lcao"???
					// mohan update 2010-09-30
					wf.LCAO_in_pw_k(ik, wf.wanf2[ik]);
                    //LOCAL_BASIS=0; xiaohui modify 2013-09-01
                    BASIS_TYPE="pw"; //xiaohui add 2013-09-01
                }
            }

            //xiaohui modify 2013-09-01. Attention! Maybe there is some problem. When will "LOCAL_BASIS==1"???
            //if (LOCAL_BASIS == 1)
            //{
            //    winput::imp_pao = 2;
            //    winput::b_recon = true;
            //    winput::sph_proj = 0;
            //    //wannier::running( ofs_running );
            //} xiaohui modify 2013-09-01

            //Spillage sp;
            //sp.get_both(NBANDS, NLOCAL, wf.wanf2, wf.evc);
        }

        // output overlap
        if ( winput::out_spillage <= 2 )
        {
            Numerical_Basis numerical_basis;
            numerical_basis.output_overlap(wf.evc);
            DONE(ofs_running,"BASIS OVERLAP (Q and S) GENERATION.");
        }
    }


	//xiaohui add 2013-09-01. Attention! Maybe there is some problem.
	if(Optical::opt_epsilon2 && (BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw"))
	{
		Optical opt;
		opt.cal_epsilon2(NBANDS);
	}

	// compute density of states
	en.perform_dos();
	timer::tick("Run_Frag","plane_wave_line",'B');

    return;
}
