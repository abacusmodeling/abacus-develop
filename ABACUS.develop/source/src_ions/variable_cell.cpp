#include "variable_cell.h"
#include "../src_pw/global.h"
#include "../input.h"
//#include "src_pw/algorithms.h"
//#include "src_pw/pseudopot_cell_us.h"
//#include "src_pw/optical.h"
//#include "src_pw/cal_test.h"
//#include "src_lcao/dftu.h"   //Quxin add for DFT+U on 20201029
//#include "src_pw/winput.h"
//#include "src_lcao/sltk_atom_arrange.h"

Variable_Cell::Variable_Cell(){}
Variable_Cell::~Variable_Cell(){}


//LiuXh add a new function here,
//which is used to do initialization after variable cell
//20180515
void Variable_Cell::init_after_vc(void)
{
	TITLE("Variable_Cell","init_after_vc");

    ucell.setup_cell_after_vc(global_pseudo_dir, global_atom_card, ofs_running);
    DONE(ofs_running, "SETUP UNITCELL");

    if(Symmetry::symm_flag)
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
    ppcell.init_vloc(pw.nggm);
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
void Variable_Cell::final_calculation_after_vc(void)
{
	TITLE("Variable_Cell","final_after_vc");

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
    if(Symmetry::symm_flag)
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
    pot.allocate(pw.nrxx);
    //=====================
    // init wave functions
    //=====================
    if(BASIS_TYPE=="pw")
    {
        wf.allocate(kv.nks);
    }
    else
    {
        wf.allocate_ekb_wg(kv.nks);
    }
    UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    if(BASIS_TYPE=="pw") ppcell.init(ucell.ntype);
    //=====================
    // init hamiltonian
    //=====================
    hm.hpw.init(wf.npwx, NPOL, ppcell.nkb, pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc(pw.nggm);
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
    pot.init_pot(0, pw.strucFac);//atomic_rho, v_of_rho, set_vrs

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

	return;
}
