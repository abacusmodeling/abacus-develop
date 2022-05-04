#include "esolver_ks_lcao.h"

//--------------temporary----------------------------
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../src_io/print_info.h"
#include "src_lcao/dftu.h"
#include "src_lcao/dmft.h"

#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif
//-----force& stress-------------------
#include "src_lcao/FORCE_STRESS.h"


//---------------------------------------------------

namespace ModuleESolver
{

ESolver_KS_LCAO::ESolver_KS_LCAO()
{
    classname = "ESolver_KS_LCAO";
    basisname = "LCAO";
}
ESolver_KS_LCAO::~ESolver_KS_LCAO()
{
    this->orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
}

void ESolver_KS_LCAO::Init(Input& inp, UnitCell_pseudo& ucell)
{    
    // setup GlobalV::NBANDS
	// Yu Liu add 2021-07-03
	GlobalC::CHR.cal_nelec();

	// it has been established that that
	// xc_func is same for all elements, therefore
	// only the first one if used
	if(ucell.atoms[0].xc_func=="HSE" || ucell.atoms[0].xc_func=="PBE0")
	{
		XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(ucell.atoms[0].xc_func);
	}

    //ucell.setup_cell( GlobalV::global_pseudo_dir , GlobalV::stru_file , GlobalV::ofs_running, GlobalV::NLOCAL, GlobalV::NBANDS);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set(GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(ucell, GlobalC::kv);

//--------------------------------------
// cell relaxation should begin here
//--------------------------------------

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << std::endl;

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
        GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04
	// Calculate Structure factor
    GlobalC::pw.setup_structure_factor();

	// Inititlize the charge density.
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT CHARGE");

	// Initializee the potential.
    GlobalC::pot.allocate(GlobalC::pw.nrxx);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT POTENTIAL");

#ifdef __MPI 
	if(GlobalV::CALCULATION=="nscf")
	{
		switch(GlobalC::exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				XC_Functional::set_xc_type(ucell.atoms[0].xc_func);
				break;
		}
	}
#endif

#ifdef __DEEPKS
	//wenfei 2021-12-19
	//if we are performing DeePKS calculations, we need to load a model
	if (GlobalV::deepks_scf)
	{
		// load the DeePKS model from deep neural network
		GlobalC::ld.load_model(INPUT.deepks_model);
	}
#endif

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    GlobalC::wf.allocate_ekb_wg(GlobalC::kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP
    GlobalC::UFFT.allocate();

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
	// without structure factors
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step=0;
    GlobalC::pot.init_pot(ion_step, GlobalC::pw.strucFac);

#ifdef __MPI  
	// PLEASE simplify the Exx_Global interface
	// mohan add 2021-03-25
	// Peize Lin 2016-12-03
	if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
	{
		switch(GlobalC::exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_lcao.init();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	}	
#endif

	// PLEASE do not use INPUT global variable
	// mohan add 2021-03-25
	// Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		GlobalC::dftu.init(ucell, this->LM);
	}

    if (INPUT.dft_plus_dmft) GlobalC::dmft.init(INPUT, ucell);
    
    //------------------init Basis_lcao----------------------
    // Init Basis should be put outside of Ensolver.
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    this->Init_Basis_lcao(this->orb_con, inp, ucell);
    //------------------init Basis_lcao----------------------

    //------------------init Hamilt_lcao----------------------
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    this->LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, orb_con.ParaV);
    //------------------init Hamilt_lcao----------------------
  //init Psi
    if (GlobalV::GAMMA_ONLY_LOCAL)
        this->LOWF.wfc_gamma.resize(GlobalV::NSPIN);
    else
        this->LOWF.wfc_k.resize(GlobalC::kv.nks);
    //pass Hamilt-pointer to Operator
    this->UHM.genH.LM = this->UHM.LM = &this->LM;
    //pass basis-pointer to EState and Psi
    this->LOC.ParaV = this->LOWF.ParaV = this->LM.ParaV;
}

void ESolver_KS_LCAO::Run(int istep, UnitCell_pseudo& cell)
{
    this->beforescf(istep);
    this->solver(istep, this->LOC, this->LOWF);
    this->afterscf(1);
    return;
}

void ESolver_KS_LCAO::cal_Energy(energy &en)
{
    
}

void ESolver_KS_LCAO::cal_Force(ModuleBase::matrix &force)
{
    
    Force_Stress_LCAO FSL(this->RA);
    FSL.getForceStress(GlobalV::CAL_FORCE, GlobalV::CAL_STRESS,
        GlobalV::TEST_FORCE, GlobalV::TEST_STRESS,
        this->LOC, this->LOWF, this->UHM, force, this->scs);
    //delete RA after cal_Force
    this->RA.delete_grid();
}
void ESolver_KS_LCAO::cal_Stress(ModuleBase::matrix &stress)
{
    //copy the stress
    stress = this->scs;
}
void ESolver_KS_LCAO::postprocess()
{
    GlobalC::en.perform_dos(this->LOWF,this->UHM);
}

void ESolver_KS_LCAO::Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell_pseudo& ucell)
{
    // Set the variables first
    this->orb_con.gamma_only = GlobalV::GAMMA_ONLY_LOCAL;
    this->orb_con.nlocal = GlobalV::NLOCAL;
    this->orb_con.nbands = GlobalV::NBANDS;
    this->orb_con.ParaV.nspin = GlobalV::NSPIN;
    this->orb_con.dsize = GlobalV::DSIZE;
    this->orb_con.nb2d = GlobalV::NB2D;
    this->orb_con.dcolor = GlobalV::DCOLOR;
    this->orb_con.drank = GlobalV::DRANK;
    this->orb_con.myrank = GlobalV::MY_RANK;
    this->orb_con.calculation = GlobalV::CALCULATION;
    this->orb_con.ks_solver = GlobalV::KS_SOLVER;
    this->orb_con.setup_2d = true;

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    this->orb_con.read_orb_first(
        GlobalV::ofs_running,
        GlobalC::ORB,
        ucell.ntype,
        ucell.lmax,
        inp.lcao_ecut,
        inp.lcao_dk,
        inp.lcao_dr,
        inp.lcao_rmax,
        GlobalV::deepks_setorb,
        inp.out_mat_r,
        GlobalV::CAL_FORCE,
        GlobalV::MY_RANK);

    ucell.infoNL.setupNonlocal(
        ucell.ntype,
		ucell.atoms,
		GlobalV::ofs_running,
        GlobalC::ORB);

#ifdef __MPI   
	this->orb_con.set_orb_tables(
		GlobalV::ofs_running,
		GlobalC::UOT,
		GlobalC::ORB,
		ucell.lat0,
		GlobalV::deepks_setorb,
		Exx_Abfs::Lmax,
		ucell.infoNL.nprojmax,
		ucell.infoNL.nproj,
        ucell.infoNL.Beta);
#else
	int Lmax=0;
	this->orb_con.set_orb_tables(
		GlobalV::ofs_running,
		GlobalC::UOT,
		GlobalC::ORB,
		ucell.lat0,
		GlobalV::deepks_setorb,
		Lmax,
		ucell.infoNL.nprojmax,
	    ucell.infoNL.nproj,
        ucell.infoNL.Beta);
#endif

    if (this->orb_con.setup_2d)
        this->orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);
}


void ESolver_KS_LCAO::eachiterinit(int iter)
{
    
}
void ESolver_KS_LCAO::hamilt2density(int istep, int iter, double ethr) 
{
    
}
void ESolver_KS_LCAO::updatepot(bool conv)
{
    
}
void ESolver_KS_LCAO::eachiterfinish(int iter, bool conv)
{
    
}
void ESolver_KS_LCAO::afterscf(bool conv_elec)
{
    //1. DeePKS calculation 
    #ifdef __DEEPKS
    const Parallel_Orbitals* pv = this->LOWF.ParaV;
    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf)
    {
        //this part is for integrated test of deepks
        //so it is printed no matter even if deepks_out_labels is not used
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            GlobalC::ld.cal_projected_DM(this->LOC.dm_gamma[0],
                GlobalC::ucell,
                GlobalC::ORB,
                GlobalC::GridD,
                pv->trace_loc_row,
                pv->trace_loc_col);
        }
        else
        {
            GlobalC::ld.cal_projected_DM_k(this->LOC.dm_k,
                GlobalC::ucell,
                GlobalC::ORB,
                GlobalC::GridD,
                pv->trace_loc_row,
                pv->trace_loc_col,
                GlobalC::kv.nks,
                GlobalC::kv.kvec_d);
        }
        GlobalC::ld.cal_descriptor();    //final descriptor
        GlobalC::ld.check_descriptor(GlobalC::ucell);
        
        if (GlobalV::deepks_out_labels) GlobalC::ld.save_npy_d(GlobalC::ucell.nat);            //libnpy needed
    }

    if (GlobalV::deepks_scf)
    {
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            GlobalC::ld.cal_e_delta_band(this->LOC.dm_gamma,
                pv->trace_loc_row,
                pv->trace_loc_col,
                pv->nrow);
        }
        else
        {
            GlobalC::ld.cal_e_delta_band_k(this->LOC.dm_k,
            pv->trace_loc_row,
            pv->trace_loc_col,
            GlobalC::kv.nks,
            pv->nrow,
            pv->ncol);
        }
        std::cout << "E_delta_band = " << std::setprecision(8) << GlobalC::ld.e_delta_band << " Ry" << " = " << std::setprecision(8) << GlobalC::ld.e_delta_band * ModuleBase::Ry_to_eV << " eV" << std::endl;
        std::cout << "E_delta_NN= "<<std::setprecision(8) << GlobalC::ld.E_delta << " Ry" << " = "<<std::setprecision(8)<<GlobalC::ld.E_delta*ModuleBase::Ry_to_eV<<" eV"<<std::endl;
    }
#endif
    //2. some outputs
    if(INPUT.dft_plus_dmft)
    {
    // Output sparse overlap matrix S(R)
    this->output_SR("outputs_to_DMFT/overlap_matrix/SR.csr");
    
    // Output wave functions, bands, k-points information, and etc.
    GlobalC::dmft.out_to_dmft(this->LOWF, *this->UHM.LM);
    }

    if(Pdiag_Double::out_mat_hsR)
    {
        this->output_HS_R(); //LiuXh add 2019-07-15
    }

}

}