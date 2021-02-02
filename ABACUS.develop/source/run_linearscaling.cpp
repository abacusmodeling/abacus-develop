#include "run_linearscaling.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_pw/algorithms.h"
#include "src_pw/pseudopot_cell_us.h"
#include "src_lcao/sltk_atom_arrange.h"
#include "src_lcao/local_orbital_ions.h"

Run_linearscaling::Run_linearscaling(){}
Run_linearscaling::~Run_linearscaling(){}


void Run_linearscaling::linear_scaling_line(void)
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
