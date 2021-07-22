#include "variable_cell.h"
#include "../src_pw/global.h"
#include "../input.h"

Variable_Cell::Variable_Cell(){}
Variable_Cell::~Variable_Cell(){}

void Variable_Cell::init_after_vc(void)
{
	TITLE("Variable_Cell","init_after_vc");

    ucell.setup_cell_after_vc(GlobalV::global_pseudo_dir, GlobalC::out, GlobalV::global_atom_card, GlobalV::ofs_running);
    DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    if(Symmetry::symm_flag)
    {
        symm.analy_sys(ucell, GlobalC::out);
        DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    GlobalC::kv.set_after_vc(symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec);
    DONE(GlobalV::ofs_running, "INIT K-POINTS");

    GlobalC::pw.update_gvectors(GlobalV::ofs_running, ucell);

    GlobalC::pw.setup_structure_factor();

    if(GlobalV::BASIS_TYPE=="pw")
    {
        wf.init_after_vc(GlobalC::kv.nks);
        wf.init_at_1();
    }

    GlobalV::ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << endl;
    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc(GlobalC::pw.nggm, ppcell.vloc);
    DONE(GlobalV::ofs_running,"LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    if(GlobalV::BASIS_TYPE=="pw")
    {
        ppcell.init_vnl(ucell);
        DONE(GlobalV::ofs_running,"NON-LOCAL POTENTIAL");
    }

    return;
}
    

void Variable_Cell::final_calculation_after_vc(void)
{
	TITLE("Variable_Cell","final_after_vc");

    cout<<" -----------------------------------------------------------------"<<endl;

    cout<<"\n -----------------------------------------------------------------"<<endl;
    cout<<" The structure has been fully relaxed or MD finished, and the following is a scf"<<endl;
    cout<<" calculation at the final structure. The fft grids and G-vectors "<<endl;
    cout<<" are recalculated for the final relaxed unit cell."<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;

    OUT(GlobalV::ofs_running," ------------------------------------------------------------------------------------");

    OUT(GlobalV::ofs_running,"\n ------------------------------------------------------------------------------------");
    OUT(GlobalV::ofs_running," The structure has been fully relaxed, and the following is a scf calculation");
    OUT(GlobalV::ofs_running," at the final structure. The fft grids and G-vectors are recalculated for the");
    OUT(GlobalV::ofs_running," final relaxed unit cell.");
    OUT(GlobalV::ofs_running," ------------------------------------------------------------------------------------");

    // (5) Setup the unitcell.
    ucell.setup_cell_after_vc(GlobalV::global_pseudo_dir, GlobalC::out, GlobalV::global_atom_card, GlobalV::ofs_running);
    DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // (6) symmetry analysize.
    if(Symmetry::symm_flag)
    {
        symm.analy_sys(ucell, GlobalC::out);
        DONE(GlobalV::ofs_running, "SYMMETRY");
    }


    // (7) Setup the k points according to symmetry.
    GlobalC::kv.set( symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec );
    DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // (1) Init the plane wave.
    GlobalC::pw.gen_pw(GlobalV::ofs_running, ucell, GlobalC::kv);
    DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << endl;
    cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << endl;

    // init the grid, then the charge
    // on grid can be distributed.
    Pgrid.init_final_scf(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, 
		GlobalC::pw.nczp, GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04

    //=====================
    // init potential
    //=====================
    CHR.init_final_scf();
    pot.allocate(GlobalC::pw.nrxx);

    //=====================
    // init wave functions
    //=====================
    if(GlobalV::BASIS_TYPE=="pw")
    {
        wf.allocate(GlobalC::kv.nks);
    }
    else
    {
        wf.allocate_ekb_wg(GlobalC::kv.nks);
    }
    GlobalC::UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    if(GlobalV::BASIS_TYPE=="pw")
	{
		ppcell.init(ucell.ntype);
	}

    //=====================
    // init hamiltonian
    //=====================
    hm.hpw.allocate(wf.npwx, GlobalV::NPOL, ppcell.nkb, GlobalC::pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc(GlobalC::pw.nggm, ppcell.vloc);
    DONE(GlobalV::ofs_running,"LOCAL POTENTIAL");
    //======================================
    // Initalize non local pseudopotential
    //======================================
    if(GlobalV::BASIS_TYPE=="pw")
    {
        ppcell.init_vnl(ucell);
        DONE(GlobalV::ofs_running,"NON-LOCAL POTENTIAL");
    }
    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    pot.init_pot(0, GlobalC::pw.strucFac);//atomic_rho, v_of_rho, set_vrs

    //==================================================
    // Create ppcell.tab_at , for trial wave functions.
    // Initial start wave functions
    //================================
    if(GlobalV::BASIS_TYPE=="pw")
    {
		pot.newd();

		wf.init_at_1();

        wf.wfcinit();
        DONE(GlobalV::ofs_running,"INIT BASIS");
    }

	return;
}
