//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-06
//==========================================================
#include "run_frag.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_pw/algorithms.h"
#include "src_pw/pseudopot_cell_us.h"
#include "src_pw/optical.h"
#include "src_pw/cal_test.h"
#include "src_lcao/dftu.h"   //Quxin add for DFT+U on 20201029
#include "src_pw/winput.h"
#include "src_lcao/sltk_atom_arrange.h"
#include "src_lcao/local_orbital_ions.h"

Run_Frag::Run_Frag(){}
Run_Frag::~Run_Frag(){}

// * initialize the NPOOL for k-points
// * initialize the input parameters for wannier functions
// * setup the unit cell
// * symmetry analysis
// * setup the k-points according to symmetry annalysis
void Run_Frag::init(void)
{
	TITLE("Run_Frag","init");
	timer::tick("Run_Frag","init",'B');

#ifdef __MPI
    // (1) If k point number > 1, After reading in NPOOL, 
	// divide the NPROC processprs into NPOOL.
    Pkpoints.init();
#endif

    // (2) Read in parameters about wannier functions.
    winput::Init( global_wannier_card );

    //xiaohui move 3 lines, 2015-09-30
    //stringstream ss2;
    //ss2 << global_out_dir << "INPUTw";
    //winput::Print( ss2.str() );

    // (3) Print the parameters into INPUT file.
    stringstream ss1;
    ss1 << global_out_dir << global_in_card;
    INPUT.Print( ss1.str() );
    //DONE(ofs_running,"READING CARDS");

    // (4) Setup the unitcell.
    ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running);
    DONE(ofs_running, "SETUP UNITCELL");

    // (5) symmetry analysize.
    if (SYMMETRY)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }
    
	// (6) Setup the k points according to symmetry.
	kv.set( symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");
	

	// (7) check the number of basis, the warning should be moved to 
	// other places -- mohan 2021-01-30
	// mohan add 2011-01-5
	if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") 
	{
		if( NLOCAL < NBANDS )
		{
			WARNING_QUIT("UnitCell_pseudo::cal_nwfc","NLOCAL < NBANDS");
		}
		else
		{
			//OUT(ofs_running,"NLOCAL",NLOCAL);
			OUT(ofs_running,"NBASE",NLOCAL);
			OUT(ofs_running,"NBANDS",NBANDS);
		}
	}

	timer::tick("Run_Frag","init",'B');
    return;
}


void Run_Frag::LCAO_line(void)
{
    TITLE("Run_Frag","lcao_line");
	timer::tick("Run_Frag","LCAO_line",'B');
	
	// (1) Inititlize the charge density.
    CHR.init();
    DONE(ofs_running,"INIT CHARGE");

	// (2) Initializee the potential.
    pot.init(pw.nrxx);
    DONE(ofs_running,"INIT POTENTIAL");

    // declration
    enum use_wf_coef {SOME_PW, ALL_LO};
    use_wf_coef uoc = ALL_LO;

	switch (uoc)
	{
		case ALL_LO:
			// (4) Init the local wave functions.
			wf.init_local();
			// (5) Init the FFT.
			UFFT.allocate();
			// (6) Init the hamiltonian. 
			// first0 stands for nkb, but no used.
			// second0 stands for no use hpw.init()
			hm.init(0);
			// (7) Init the local part of NC pseudopotential.
			ppcell.init_vloc();
			// (8) Init the potential.
			pot.init_pot(0);//atomic_rho, v_of_rho, set_vrs
			break;
		case SOME_PW:
			wf.init(kv.nks);
			UFFT.allocate();
			ppcell.init(ucell.ntype);
			hm.init();
			ppcell.init_vloc();
			ppcell.init_vnl();
			pot.init_pot(0);//atomic_rho, v_of_rho, set_vrs
			pot.newd();//once
			DONE(ofs_running,"INIT POTENTIAL");
			wf.wfcinit();
			DONE(ofs_running,"INIT SOME_PW");
			break;
	}

    // Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		dftu.init();
	}

	Local_Orbital_Ions ions;
	ions.opt_ions();
	en.perform_dos();

	timer::tick("Run_Frag","LCAO_line",'B');
    return;
}



void Run_Frag::plane_wave_line(void)
{
    TITLE("Run_Frag","plane_wave_line");
	timer::tick("Run_Frag","plane_wave_line",'B');

//----------------------------------------------------------
// 1 read in initial data:
//   a lattice structure:atom_species,atom_positions,lattice vector
//   b k_points
//   c pseudopotential
// 2 setup planeware basis, FFT,structure factor, ...
// 3 initialize local and nonlocal pseudopotential in G_space
// 4 initialize charge desity and warefunctios in G_space
//----------------------------------------------------------

    //=====================
    // init potential
    //=====================
    CHR.init();
    pot.init(pw.nrxx);

    //=====================
    // init wave functions
    //=====================
    wf.init(kv.nks);
	UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    ppcell.init(ucell.ntype);

    //=====================
    // init hamiltonian
    //=====================
    hm.init();
//  DONE(ofs_running,"CHARGE, POTENTIAL, WAVE FUNCTINOS ALLOCATION");

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
   	wf.wfcinit();


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


void Run_Frag::linear_scaling_line(void)
{
	TITLE("Run_Frag","linear_scaling_line");

    // (2) Init the charge density.
    CHR.init();
    DONE(ofs_running,"INIT CHARGE");

    // (3) Init the potential.
    pot.init(pw.nrxx);
    DONE(ofs_running,"INIT POTENTIAL");

    // declration
    enum use_wf_coef {SOME_PW, ALL_LO};
    // generate object
    use_wf_coef uoc = ALL_LO;

    switch (uoc)
    {
    case ALL_LO:

        // (4) Init the local wave functions.
        wf.init_local();

        // (5) Init the FFT.
        UFFT.allocate();

        // (6) Init the hamiltonian.
        hm.init();

        // (7) Init the local part of NC pseudopotential.
        ppcell.init_vloc();

        // (8) Init the potential.
        pot.init_pot(0);//atomic_rho, v_of_rho, set_vrs
        DONE(ofs_running,"INIT ALL_LO");
        break;

		case SOME_PW:
        wf.init(kv.nks);
        UFFT.allocate();
        ppcell.init(ucell.ntype);
        hm.init();
        ppcell.init_vloc();
        ppcell.init_vnl();
        pot.init_pot(0);//atomic_rho, v_of_rho, set_vrs
        pot.newd();//once
        DONE(ofs_running,"INIT POTENTIAL");
        wf.wfcinit();
        DONE(ofs_running,"INIT SOME_PW");
        break;
    }

    // (9) Begin the ion iteration.
    Local_Orbital_Ions ions;
   	ions.opt_ions();

    return;
}


//LiuXh add a new function here,
//which is used to do initialization after variable cell
//20180515
void Run_Frag::init_after_vc(void)
{
	TITLE("Run_Frag","init_after_vc");

    ucell.setup_cell_after_vc(global_pseudo_dir, global_atom_card, ofs_running);
    DONE(ofs_running, "SETUP UNITCELL");

    if(SYMMETRY)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }

    kv.set_after_vc(symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec);
    DONE(ofs_running, "INIT K-POINTS");

    pw.update_gvectors(ofs_running, ucell);

    pw.setup_structure_factor();

    if(BASIS_TYPE=="pw")
    {
        wf.init_after_vc(kv.nks);
        wf.init_at_1();
    }

    ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << endl;
    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc();
    DONE(ofs_running,"LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    if(BASIS_TYPE=="pw")
    {
        ppcell.init_vnl();
        DONE(ofs_running,"NON-LOCAL POTENTIAL");
    }

/*
    pot.init_pot(0);

    ofs_running << " Setup the new wave functions?" << endl;
    wf.wfcinit();
*/

    return;
}
    

//LiuXh add a new function here,
//which is used to do initialization after variable cell
//20180619
void Run_Frag::final_calculation_after_vc(void)
{
    //cout<<" -------------------------------------------"<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;

    cout<<"\n -----------------------------------------------------------------"<<endl;
    cout<<" The structure has been fully relaxed, and the following is a scf"<<endl;
    cout<<" calculation at the final structure. The fft grids and G-vectors "<<endl;
    cout<<" are recalculated for the final relaxed unit cell."<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;

    OUT(ofs_running," ------------------------------------------------------------------------------------");

    OUT(ofs_running,"\n ------------------------------------------------------------------------------------");
    OUT(ofs_running," The structure has been fully relaxed, and the following is a scf calculation");
    OUT(ofs_running," at the final structure. The fft grids and G-vectors are recalculated for the");
    OUT(ofs_running," final relaxed unit cell.");
    OUT(ofs_running," ------------------------------------------------------------------------------------");

    // (5) Setup the unitcell.
    ucell.setup_cell_after_vc(global_pseudo_dir, global_atom_card, ofs_running);
    DONE(ofs_running, "SETUP UNITCELL");

    // (6) symmetry analysize.
    if (SYMMETRY)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }


    // (7) Setup the k points according to symmetry.
    kv.set( symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");

    // (1) Init the plane wave.
    pw.gen_pw(ofs_running, ucell, kv);
    DONE(ofs_running,"INIT PLANEWAVE");
    cout << " UNIFORM GRID DIM     : " << pw.nx <<" * " << pw.ny <<" * "<< pw.nz << endl;
    cout << " UNIFORM GRID DIM(BIG): " << pw.nbx <<" * " << pw.nby <<" * "<< pw.nbz << endl;

    // init the grid, then the charge
    // on grid can be distributed.
    Pgrid.init_final_scf(pw.ncx, pw.ncy, pw.ncz, pw.nczp, pw.nrxx, pw.nbz, pw.bz); // mohan add 2010-07-22, update 2011-05-04

    //=====================
    // init potential
    //=====================
    CHR.init_final_scf();
    pot.init(pw.nrxx);
    //=====================
    // init wave functions
    //=====================
    if(BASIS_TYPE=="pw")
    {
        wf.init(kv.nks);
    }
    else
    {
        wf.init_local();
    }
    UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    if(BASIS_TYPE=="pw") ppcell.init(ucell.ntype);
    //=====================
    // init hamiltonian
    //=====================
    hm.init();

    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc();
    DONE(ofs_running,"LOCAL POTENTIAL");
    //======================================
    // Initalize non local pseudopotential
    //======================================
    if(BASIS_TYPE=="pw")
    {
        ppcell.init_vnl();
        DONE(ofs_running,"NON-LOCAL POTENTIAL");
    }
    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    pot.init_pot(0);//atomic_rho, v_of_rho, set_vrs

    if(BASIS_TYPE=="pw") pot.newd();//once
    DONE(ofs_running,"INIT POTENTIAL");

    //==================================================
    // create ppcell.tab_at , for trial wave functions.
    //==================================================
    if(BASIS_TYPE=="pw") wf.init_at_1();
    //================================
    // Initial start wave functions
    //================================
    if(BASIS_TYPE=="pw")
    {
        wf.wfcinit();
        DONE(ofs_running,"INIT BASIS");
    }
}
