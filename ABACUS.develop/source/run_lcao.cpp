//==========================================================
// AUTHOR : mohan
// DATE : 2021-02-01
//==========================================================
#include "run_lcao.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_lcao/dftu.h"   //Quxin add for DFT+U on 20201029
#include "src_io/winput.h"
#include "src_global/sltk_atom_arrange.h"
#include "src_lcao/local_orbital_ions.h"

Run_lcao::Run_lcao(){}
Run_lcao::~Run_lcao(){}


void Run_lcao::lcao_line(void)
{
    TITLE("Run_lcao","lcao_line");
	timer::tick("Run_lcao","lcao_line",'B');
	
	// (1) Inititlize the charge density.
    CHR.init();
    DONE(ofs_running,"INIT CHARGE");

	// (2) Initializee the potential.
    pot.init(pw.nrxx);
    DONE(ofs_running,"INIT POTENTIAL");

    // declration
    enum use_wf_coef {SOME_PW, ALL_LO};
    use_wf_coef uoc = ALL_LO;

	// Peize Lin add 2018-11-30
	if(CALCULATION=="nscf")
	{
		switch(exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_global.info.set_xcfunc(xcf);
				break;
		}
	}

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

    // Peize Lin 2016-12-03
	if (CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax" || CALCULATION=="cell-relax")
	{
		switch(exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_lcao.init();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		}
	}	

    // Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		dftu.init();
	}

	Local_Orbital_Ions ions;
	ions.opt_ions();
	en.perform_dos();

	timer::tick("Run_lcao","lcao_line",'B');
    return;
}
