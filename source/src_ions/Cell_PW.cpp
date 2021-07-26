#include "Cell_PW.h"
#include "ions.h"

void Cell_PW::opt_cells_pw()
{
    TITLE("Cell_PW", "opt_cells_pw");
    timer::tick("Cell_PW", "opt_cells_pw");
    GlobalC::wf.allocate(GlobalC::kv.nks);

    GlobalC::UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    GlobalC::ppcell.init(GlobalC::ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);
    DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    pot.init_pot(0, GlobalC::pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    pot.newd();

    DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    //==================================================
    // create GlobalC::ppcell.tab_at , for trial wave functions.
    //==================================================
    GlobalC::wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if (GlobalV::NBANDS != 0 || (GlobalV::CALCULATION != "scf-sto" && GlobalV::CALCULATION != "relax-sto" && GlobalV::CALCULATION != "md-sto")) //qianrui add
    {
        GlobalC::wf.wfcinit();
    }
#ifdef __LCAO
    switch (GlobalC::exx_global.info.hybrid_type) // Peize Lin add 2019-03-09
    {
    case Exx_Global::Hybrid_Type::HF:
    case Exx_Global::Hybrid_Type::PBE0:
    case Exx_Global::Hybrid_Type::HSE:
        GlobalC::exx_lip.init(&GlobalC::kv, &GlobalC::wf, &GlobalC::pw, &GlobalC::UFFT, &GlobalC::ucell);
        break;
    case Exx_Global::Hybrid_Type::No:
        break;
    case Exx_Global::Hybrid_Type::Generate_Matrix:
    default:
        throw invalid_argument(TO_STRING(__FILE__) + TO_STRING(__LINE__));
    }
#endif

    DONE(GlobalV::ofs_running, "INIT BASIS");

    // ion optimization begins
    // electron density optimization is included in ion optimization

    Ions ions;
    ions.opt_ions_pw();
    
    timer::tick("Cell_PW", "opt_cells_pw");
}
