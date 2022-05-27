#include "variable_cell.h"
#include "../src_pw/global.h"
#include "../input.h"

Variable_Cell::Variable_Cell(){}
Variable_Cell::~Variable_Cell(){}

void Variable_Cell::init_after_vc(ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("Variable_Cell","init_after_vc");

    GlobalC::ucell.setup_cell_after_vc(GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    if(ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    GlobalC::kv.set_after_vc(GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    GlobalC::pw.update_gvectors(GlobalV::ofs_running, GlobalC::ucell);

    GlobalC::pw.setup_structure_factor();

    if(GlobalV::BASIS_TYPE=="pw")
    {
        GlobalC::wf.init_after_vc(GlobalC::kv.nks, p_esolver->psi);
        GlobalC::wf.init_at_1();
    }

    GlobalV::ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << std::endl;
    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    if(GlobalV::BASIS_TYPE=="pw")
    {
        GlobalC::ppcell.init_vnl(GlobalC::ucell);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"NON-LOCAL POTENTIAL");
    }

    return;
}
