#include "esolver_ks_lcao.h"
#include "module_io/cal_r_overlap_R.h"
#include "module_io/dm_io.h"
#include "module_io/write_dm_sparse.h"
#include "module_io/rho_io.h"
#include "module_io/write_HS_R.h"
#include "module_io/write_dos_lcao.h"
#include "module_io/write_istate_info.h"
#include "module_io/mulliken_charge.h"
#include "module_io/nscf_band.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_io/nscf_fermi_surf.h"

//--------------temporary----------------------------
#include "module_base/global_function.h"
#include "module_io/print_info.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_elecstate/occupy.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#ifdef __EXX
// #include "module_rpa/rpa.h"
#include "module_ri/RPA_LRI.h"
#endif

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
// function used by deepks
#include "module_elecstate/cal_dm.h"
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

void ESolver_KS_LCAO::Init(Input& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "Init");
    // if we are only calculating S, then there is no need
    // to prepare for potentials and so on

    if (GlobalV::CALCULATION == "get_S")
    {
        ucell.read_pseudo(GlobalV::ofs_running);

        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            GlobalC::symm.analy_sys(ucell, GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
        }

        // Setup the k points according to symmetry.
        GlobalC::kv.set(GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

        Print_Info::setup_parameters(ucell, GlobalC::kv);
    }
    else
    {
        ESolver_KS::Init(inp, ucell);
    } // end ifnot get_S

    //init ElecState
    //autoset nbands in ElecState, it should before basis_init (for Psi 2d divid)
    if(this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecStateLCAO(&(chr),
                                                   &(GlobalC::kv),
                                                   GlobalC::kv.nks,
                                                   &(this->LOC),
                                                   &(this->UHM),
                                                   &(this->LOWF));
    }

    //------------------init Basis_lcao----------------------
    // Init Basis should be put outside of Ensolver.
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    this->Init_Basis_lcao(this->orb_con, inp, ucell);
    //------------------init Basis_lcao----------------------

    // pass Hamilt-pointer to Operator
    this->UHM.genH.LM = this->UHM.LM = &this->LM;
    // pass basis-pointer to EState and Psi
    this->LOC.ParaV = this->LOWF.ParaV = this->LM.ParaV = &(this->orb_con.ParaV);

    if (GlobalV::CALCULATION == "get_S")
    {
        return;
    }

    //------------------init Hamilt_lcao----------------------
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    this->LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, orb_con.ParaV);
    //------------------init Hamilt_lcao----------------------

#ifdef __EXX
    // PLEASE simplify the Exx_Global interface
    // mohan add 2021-03-25
    // Peize Lin 2016-12-03
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax"
        || GlobalV::CALCULATION == "cell-relax" || GlobalV::CALCULATION == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            /* In the special "two-level" calculation case,
            first scf iteration only calculate the functional without exact exchange.
            but in "nscf" calculation, there is no need of "two-level" method. */
            if (ucell.atoms[0].ncpp.xc_func == "HSE" || ucell.atoms[0].ncpp.xc_func == "PBE0")
            {
                XC_Functional::set_xc_type("pbe");
            }
            else if (ucell.atoms[0].ncpp.xc_func == "SCAN0")
            {
                XC_Functional::set_xc_type("scan");
            }

			// GlobalC::exx_lcao.init();
            if(GlobalC::exx_info.info_ri.real_number)
                GlobalC::exx_lri_double.init(MPI_COMM_WORLD);
            else
                GlobalC::exx_lri_complex.init(MPI_COMM_WORLD);
        }
    }
#endif

    // Quxin added for DFT+U
    if (GlobalV::dft_plus_u)
    {
        GlobalC::dftu.init(ucell, this->LM);
    }

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
    // without structure factors
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, GlobalC::rhopw);

    // init HSolver
    if(this->phsol == nullptr)
    {
        this->phsol = new hsolver::HSolverLCAO(this->LOWF.ParaV);
        this->phsol->method = GlobalV::KS_SOLVER;
    }

    // Inititlize the charge density.
    this->pelec->charge->allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);

        // Initialize the potential.
    if(this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(
            GlobalC::rhopw,
            &GlobalC::ucell,
            &(GlobalC::ppcell.vloc),
            &(GlobalC::sf.strucFac),
            &(GlobalC::en.etxc),
            &(GlobalC::en.vtxc)
        );
    }

#ifdef __DEEPKS
    // wenfei 2021-12-19
    // if we are performing DeePKS calculations, we need to load a model
    if (GlobalV::deepks_scf)
    {
        // load the DeePKS model from deep neural network
        GlobalC::ld.load_model(INPUT.deepks_model);
    }
#endif

    //Fix pelec->wg by ocp_kb
    if(GlobalV::ocp)
    {
        this->pelec->fixed_weights(GlobalV::ocp_kb.data());
    }
}

void ESolver_KS_LCAO::cal_Energy(double& etot)
{
    etot = GlobalC::en.etot;
}

void ESolver_KS_LCAO::cal_Force(ModuleBase::matrix& force)
{
    Force_Stress_LCAO FSL(this->RA);
    FSL.getForceStress(GlobalV::CAL_FORCE,
                       GlobalV::CAL_STRESS,
                       GlobalV::TEST_FORCE,
                       GlobalV::TEST_STRESS,
                       this->LOC,
                       this->pelec,
                       this->psid,
                       this->psi,
                       this->UHM,
                       force,
                       this->scs);
    // delete RA after cal_Force
    this->RA.delete_grid();
    this->have_force = true;
}

void ESolver_KS_LCAO::cal_Stress(ModuleBase::matrix& stress)
{
    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_Force(fcs);
    }
    stress = this->scs; // copy the stress
    this->have_force = false;
}

void ESolver_KS_LCAO::postprocess()
{
    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    if (GlobalC::en.out_dos != 0 || GlobalC::en.out_band != 0 || GlobalC::en.out_proj_band != 0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << std::endl;
        GlobalV::ofs_running << " | done here.                                                         |" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }
    // qianrui modify 2020-10-18
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax")
    {
        ModuleIO::write_istate_info(this->pelec,&(GlobalC::kv),&(GlobalC::Pkpoints));
    }

    // GlobalV::mulliken charge analysis
#ifdef __LCAO
    if (GlobalV::out_mul == 1)
    {
        Mulliken_Charge MC(psid, psi);
        MC.stdout_mulliken(this->UHM, this->pelec->wg);
    } // qifeng add 2019/9/10
#endif

    int nspin0 = 1;
    if (GlobalV::NSPIN == 2) nspin0 = 2;

    if (GlobalC::en.out_band) // pengfei 2014-10-13
    {
        int nks = 0;
        if (nspin0 == 1)
        {
            nks = GlobalC::kv.nkstot;
        }
        else if (nspin0 == 2)
        {
            nks = GlobalC::kv.nkstot / 2;
        }

        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << GlobalV::global_out_dir << "BANDS_" << is + 1 << ".dat";
            GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
            ModuleIO::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, GlobalC::en.ef*0, this->pelec->ekb,&(GlobalC::kv),&(GlobalC::Pkpoints));
        }
    } // out_band

    if (GlobalC::en.out_proj_band) // Projeced band structure added by jiyy-2022-4-20
    {
        ModuleIO::write_proj_band_lcao(this->psid,this->psi,this->UHM,this->pelec,&(GlobalC::kv),GlobalC::ucell,GlobalC::ORB,GlobalC::GridD);
    }

    if (GlobalC::en.out_dos)
    {
        ModuleIO::write_dos_lcao(this->psid,
            this->psi,
            this->UHM,
            this->pelec->ekb,
            this->pelec->wg,
            GlobalC::en.dos_edelta_ev,
            GlobalC::en.dos_scale,
            GlobalC::en.bcoeff);

        if (GlobalC::en.out_dos == 3)
        {
            for (int i = 0; i < nspin0; i++)
            {
                std::stringstream ss3;
                ss3 << GlobalV::global_out_dir << "Fermi_Surface_" << i << ".bxsf";
                ModuleIO::nscf_fermi_surface(ss3.str(),
                    GlobalC::kv.nks,
                    GlobalV::NBANDS,
                    GlobalC::en.ef,
                    &(GlobalC::kv),
                    &(GlobalC::Pkpoints),
                    &(GlobalC::ucell),
                    this->pelec->ekb);
            }
        }
        
        if (nspin0 == 1)
        {
            GlobalV::ofs_running << " Fermi energy is " << GlobalC::en.ef << " Rydberg" << std::endl;
        }
        else if (nspin0 == 2)
        {
            GlobalV::ofs_running << " Fermi energy (spin = 1) is " << GlobalC::en.ef_up << " Rydberg" << std::endl;
            GlobalV::ofs_running << " Fermi energy (spin = 2) is " << GlobalC::en.ef_dw << " Rydberg" << std::endl;
        }
    }
}

void ESolver_KS_LCAO::Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell& ucell)
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
    this->orb_con.read_orb_first(GlobalV::ofs_running,
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

    ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, GlobalC::ORB);

    int Lmax = 0;
#ifdef __EXX
    Lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif
    this->orb_con.set_orb_tables(GlobalV::ofs_running,
                                 GlobalC::UOT,
                                 GlobalC::ORB,
                                 ucell.lat0,
                                 GlobalV::deepks_setorb,
                                 Lmax,
                                 ucell.infoNL.nprojmax,
                                 ucell.infoNL.nproj,
                                 ucell.infoNL.Beta);

    if (this->orb_con.setup_2d)
        this->orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);
}

void ESolver_KS_LCAO::eachiterinit(const int istep, const int iter)
{

    // mohan add 2010-07-16
    // used for pulay mixing.
    if (iter == 1) GlobalC::CHR_MIX.reset();

    // mohan update 2012-06-05
    GlobalC::en.deband_harris = GlobalC::en.delta_e(this->pelec);

    // mohan move it outside 2011-01-13
    // first need to calculate the weight according to
    // electrons number.

    if (GlobalC::wf.init_wfc == "file")
    {
        if (iter == 1)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.

            if (this->psi != nullptr)
            {
                this->pelec->psiToRho(this->psi[0]);
            }
            else
            {
                this->pelec->psiToRho(this->psid[0]);
            }

            // calculate the local potential(rho) again.
            // the grid integration will do in later grid integration.

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // a puzzle remains here.
            // if I don't renew potential,
            // The scf_thr is very small.
            // OneElectron, Hartree and
            // Exc energy are all correct
            // except the band energy.
            //
            // solved by mohan 2010-09-10
            // there are there rho here:
            // rho1: formed by read in orbitals.
            // rho2: atomic rho, used to construct H
            // rho3: generated by after diagonalize
            // here converged because rho3 and rho1
            // are very close.
            // so be careful here, make sure
            // rho1 and rho2 are the same rho.
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
            this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
            GlobalC::en.delta_escf(this->pelec);
        }
    }

#ifdef __EXX
    // calculate exact-exchange
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (!GlobalC::exx_info.info_global.separate_loop && this->two_level_step)
        {
            //GlobalC::exx_lcao.cal_exx_elec(this->LOC, this->LOWF.wfc_k_grid);
            if(GlobalC::exx_info.info_ri.real_number)
                GlobalC::exx_lri_double.cal_exx_elec(this->LOC, *this->LOWF.ParaV);
            else
                GlobalC::exx_lri_complex.cal_exx_elec(this->LOC, *this->LOWF.ParaV);
        }
    }
#endif

    if (GlobalV::dft_plus_u)
    {
        GlobalC::dftu.cal_slater_UJ(pelec->charge->rho); // Calculate U and J if Yukawa potential is used
    }

#ifdef __DEEPKS
    // the density matrixes of DeePKS have been updated in each iter 
    GlobalC::ld.set_hr_cal(true);
#endif

    if(!GlobalV::GAMMA_ONLY_LOCAL)
    {
        if(this->UHM.GK.get_spin() != -1)
        {
            int start_spin = -1;
            this->UHM.GK.reset_spin(start_spin);
            this->UHM.GK.destroy_pvpR();
            this->UHM.GK.allocate_pvpR();
        }
    }
}
void ESolver_KS_LCAO::hamilt2density(int istep, int iter, double ethr)
{
    // save input rho
    pelec->charge->save_rho_before_sum_band();

    // using HSolverLCAO::solve()
    if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec->eband = 0.0;
        this->pelec->demet = 0.0;
        this->pelec->ef = 0.0;
        GlobalC::en.ef_up = 0.0;
        GlobalC::en.ef_dw = 0.0;
        if (this->psi != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, GlobalV::KS_SOLVER);
        }
        else if (this->psid != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psid[0], this->pelec, GlobalV::KS_SOLVER);
        }

        // transform energy for print
        GlobalC::en.eband = this->pelec->eband;
        GlobalC::en.demet = this->pelec->demet;
        GlobalC::en.ef = this->pelec->ef;
        if (GlobalV::out_bandgap)
        {
            if (!GlobalV::TWO_EFERMI)
            {
                GlobalC::en.cal_bandgap(this->pelec);
            }
            else
            {
                GlobalC::en.cal_bandgap_updw(this->pelec);
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
    }

    // print ekb for each k point and each band
    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        this->pelec->print_band(ik, GlobalC::en.printe, iter);
    }

#ifdef __EXX
    // Peize Lin add 2020.04.04
    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        // add exx
        // Peize Lin add 2016-12-03
        GlobalC::en.set_exx();

        if (GlobalC::restart.info_load.load_H && GlobalC::restart.info_load.load_H_finish
            && !GlobalC::restart.info_load.restart_exx)
        {
            XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
            //GlobalC::exx_lcao.cal_exx_elec(this->LOC, this->LOWF.wfc_k_grid);
            if(GlobalC::exx_info.info_ri.real_number)
                GlobalC::exx_lri_double.cal_exx_elec(this->LOC, *this->LOWF.ParaV);
            else
                GlobalC::exx_lri_complex.cal_exx_elec(this->LOC, *this->LOWF.ParaV);
            GlobalC::restart.info_load.restart_exx = true;
        }
    }
    else
    {
        GlobalC::en.exx = 0.;
    }
#endif

    // if DFT+U calculation is needed, this function will calculate
    // the local occupation number matrix and energy correction
    if (GlobalV::dft_plus_u)
    {
        if(GlobalC::dftu.omc!=2)
        {
            if (GlobalV::GAMMA_ONLY_LOCAL)
                GlobalC::dftu.cal_occup_m_gamma(iter, this->LOC.dm_gamma);
            else
                GlobalC::dftu.cal_occup_m_k(iter, this->LOC.dm_k);
        }
        GlobalC::dftu.cal_energy_correction(istep);
        GlobalC::dftu.output();
    }

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        const Parallel_Orbitals* pv = this->LOWF.ParaV;
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            GlobalC::ld.cal_e_delta_band(this->LOC.dm_gamma, pv->trace_loc_row, pv->trace_loc_col, pv->nrow);
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
    }
#endif
    // (4) mohan add 2010-06-24
    // using new charge density.
    GlobalC::en.calculate_harris();

    // (5) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(pelec->charge), GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
    }

    // (6) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization(pelec->charge, pelec->nelec_spin.data());

    // (7) calculate delta energy
    GlobalC::en.deband = GlobalC::en.delta_e(this->pelec);
}
void ESolver_KS_LCAO::updatepot(const int istep, const int iter)
{
    // (9) Calculate new potential according to new Charge Density.

    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (GlobalV::out_pot < 0) // mohan add 2011-10-10
        {
            GlobalV::out_pot = -2;
        }
    }
    if (!this->conv_elec)
    {
        if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
        this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
        GlobalC::en.delta_escf(this->pelec);
    }
    else
    {
        GlobalC::en.cal_converged(this->pelec);
    }
}
void ESolver_KS_LCAO::eachiterfinish(int iter)
{
    //-----------------------------------
    // save charge density
    //-----------------------------------
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.save_disk(*this->UHM.LM, "charge", is, pelec->charge->rho);
        }
    }
    //-----------------------------------
    // output charge density for tmp
    //-----------------------------------
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        const int precision = 3;

        std::stringstream ssc;
        std::stringstream ss1;
        ssc << GlobalV::global_out_dir << "tmp"
            << "_SPIN" << is + 1 << "_CHG";
        ModuleIO::write_rho(pelec->charge->rho_save[is], is, iter, ssc.str(), precision); // mohan add 2007-10-17

        std::stringstream ssd;

        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ssd << GlobalV::global_out_dir << "tmp"
                << "_SPIN" << is + 1 << "_DM";
        }
        else
        {
            ssd << GlobalV::global_out_dir << "tmp"
                << "_SPIN" << is + 1 << "_DM_R";
        }
        ModuleIO::write_dm(is, iter, ssd.str(), precision, this->LOC.out_dm, this->LOC.DM);
    }

    if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        const int precision = 3;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            std::stringstream ss1;
            ssc << GlobalV::global_out_dir << "tmp"
                << "_SPIN" << is + 1 << "_TAU";
            ModuleIO::write_rho(pelec->charge->kin_r_save[is], is, iter, ssc.str(), precision); // mohan add 2007-10-17
        }
    }

    // (11) calculate the total energy.
    GlobalC::en.calculate_etot();
}

void ESolver_KS_LCAO::afterscf(const int istep)
{
    // Temporary liuyu add 2022-11-07
    CE.update_all_pos(GlobalC::ucell);

    // if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    // {
    //--------------------------------------
    // 1. output charge density for converged,
    // 0 means don't need to consider iter,
    //--------------------------------------

    double** dm2d;
    if(this->LOC.out_dm1 == 1)
    {
        dm2d = new double*[GlobalV::NSPIN];
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            dm2d[is] = new double[this->LOC.ParaV->nnr];
            ModuleBase::GlobalFunc::ZEROS(dm2d[is], this->LOC.ParaV->nnr);
        }
        this->LOC.cal_dm_R(this->LOC.dm_k,this->RA,dm2d);
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        const int precision = 3;

        std::stringstream ssc;
        std::stringstream ss1;
        ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG";
        ModuleIO::write_rho(pelec->charge->rho_save[is], is, 0, ssc.str()); // mohan add 2007-10-17

        std::stringstream ssd;
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
        }
        else
        {
            ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM_R";
        }
        ModuleIO::write_dm(is, 0, ssd.str(), precision, this->LOC.out_dm, this->LOC.DM);
        if(this->LOC.out_dm1 == 1)
        {
            ModuleIO::write_dm1(is, istep, dm2d, this->LOC.ParaV, this->LOC.DMR_sparse);
        }
/* Broken, please fix it
        if (GlobalV::out_pot == 1) // LiuXh add 20200701
        {
            std::stringstream ssp;
            ssp << GlobalV::global_out_dir << "SPIN" << is + 1 << "_POT";
            this->pelec->pot->write_potential(is, 0, ssp.str(), this->pelec->pot->get_effective_v(), precision);
        }
*/
    }

    if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            std::stringstream ss1;
            ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_TAU";
            ModuleIO::write_rho(pelec->charge->kin_r_save[is], is, 0, ssc.str()); // mohan add 2007-10-17
        }
    }

    if(this->LOC.out_dm1 == 1)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] dm2d[is];
        }
        delete[] dm2d;
    }

#ifdef __EXX
    if (GlobalC::exx_info.info_global.cal_exx)                         // Peize Lin add if 2022.11.14
    {
        const std::string file_name_exx = GlobalV::global_out_dir + "HexxR_" + std::to_string(GlobalV::MY_RANK);
        if(GlobalC::exx_info.info_ri.real_number)
            GlobalC::exx_lri_double.write_Hexxs(file_name_exx);
        else
            GlobalC::exx_lri_complex.write_Hexxs(file_name_exx);
    }
#endif
    if (GlobalV::out_pot == 2)
    {
        std::stringstream ssp;
        std::stringstream ssp_ave;
        ssp << GlobalV::global_out_dir << "ElecStaticPot";
        ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
        this->pelec->pot->write_elecstat_pot(ssp.str(), ssp_ave.str(), GlobalC::rhopw, pelec->charge); //output 'Hartree + local pseudopot'
    }

    if (this->conv_elec)
    {
        GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
        GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }

    if (GlobalV::OUT_LEVEL != "m")
    {
        this->pelec->print_eigenvalue(GlobalV::ofs_running);
    }

    if (this->conv_elec)
    {
        // xiaohui add "OUT_LEVEL", 2015-09-16
        if (GlobalV::OUT_LEVEL != "m")
            GlobalV::ofs_running << std::setprecision(16);
        if (GlobalV::OUT_LEVEL != "m")
            GlobalV::ofs_running << " EFERMI = " << GlobalC::en.ef * ModuleBase::Ry_to_eV << " eV" << std::endl;
        if (GlobalV::OUT_LEVEL == "ie")
        {
            GlobalV::ofs_running << " " << GlobalV::global_out_dir << " final etot is "
                                 << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
        }
    }
    else
    {
        GlobalV::ofs_running << " !! convergence has not been achieved @_@" << std::endl;
        if (GlobalV::OUT_LEVEL == "ie" || GlobalV::OUT_LEVEL == "m") // xiaohui add "m" option, 2015-09-16
            std::cout << " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!" << std::endl;
    }

#ifdef __DEEPKS
    // calculating deepks correction to bandgap
    // and save the results

    if (GlobalV::deepks_out_labels) // caoyu add 2021-06-04
    {
        GlobalC::ld.save_npy_e(GlobalC::en.etot, "e_tot.npy");
        if (GlobalV::deepks_scf)
        {
            GlobalC::ld.save_npy_e(GlobalC::en.etot - GlobalC::ld.E_delta,
                                   "e_base.npy"); // ebase :no deepks E_delta including
        }
        else // deepks_scf = 0; base calculation
        {
            GlobalC::ld.save_npy_e(GlobalC::en.etot, "e_base.npy"); // no scf, e_tot=e_base
        }

        if (GlobalV::deepks_bandgap)
        {
            int nocc = GlobalV::nelec / 2;
            int nks = GlobalC::kv.nks;
            ModuleBase::matrix deepks_bands;
            deepks_bands.create(nks,1);
            for (int iks=0; iks<nks; iks++)
            {
                for (int hl=0; hl < 1; hl++)
                {
                    deepks_bands(iks,hl) = this->pelec->ekb(iks, nocc+hl) - this->pelec->ekb(iks, nocc-1+hl);
                }   
            }
            GlobalC::ld.save_npy_o(deepks_bands, "o_tot.npy", nks);
            if (GlobalV::deepks_scf)
            {
                int nocc = GlobalV::nelec / 2;
                ModuleBase::matrix wg_hl;
                if (GlobalV::GAMMA_ONLY_LOCAL)
                {
                    wg_hl.create(GlobalV::NSPIN, GlobalV::NBANDS);
                    std::vector<std::vector<ModuleBase::matrix>> dm_bandgap_gamma;
                    dm_bandgap_gamma.resize(GlobalV::NSPIN);
                    for (int is = 0; is < GlobalV::NSPIN; is++)
                    {
                        for (int ib = 0; ib < 1; ib++)
                        {
                            wg_hl.zero_out();
                            wg_hl(is, ib+nocc-1) = -1.0;
                            wg_hl(is, ib+nocc) = 1.0;
                            dm_bandgap_gamma[ib].resize(GlobalV::NSPIN);
                            elecstate::cal_dm(this->LOWF.ParaV, wg_hl, this->psid[0], dm_bandgap_gamma[ib]);
                        }
                    }

                    GlobalC::ld.cal_orbital_precalc(dm_bandgap_gamma,
                                                    GlobalC::ucell.nat,
                                                    GlobalC::ucell,
                                                    GlobalC::ORB,
                                                    GlobalC::GridD,
                                                    *this->LOWF.ParaV);
                    

                    GlobalC::ld.save_npy_orbital_precalc(GlobalC::ucell.nat, nks);
                    GlobalC::ld.cal_o_delta(dm_bandgap_gamma, *this->LOWF.ParaV);                    
                    GlobalC::ld.save_npy_o(deepks_bands - GlobalC::ld.o_delta,"o_base.npy", nks);
                } // end deepks_scf gamma-only;
                else // multi-k bandgap label
                {
                    wg_hl.create(GlobalC::kv.nks, GlobalV::NBANDS);
                    std::vector<std::vector<ModuleBase::ComplexMatrix>> dm_bandgap_k;
                    dm_bandgap_k.resize(1);

                    for (int ib = 0; ib < 1; ib++)
                    {
                        wg_hl.zero_out();
                        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
                        {
                            wg_hl(ik, ib+nocc-1) = -1.0;
                            wg_hl(ik, ib+nocc) = 1.0;
                        }
                        dm_bandgap_k[ib].resize(GlobalC::kv.nks);
                        elecstate::cal_dm(this->LOWF.ParaV, wg_hl, this->psi[0], dm_bandgap_k[ib]);
                    }
                    
                    //GlobalC::ld.cal_o_delta_k(dm_bandgap_k, *this->LOWF.ParaV, GlobalC::kv.nks);
                    GlobalC::ld.cal_orbital_precalc_k(dm_bandgap_k,
                                                      GlobalC::ucell.nat,
                                                      GlobalC::kv.nks,
                                                      GlobalC::kv.kvec_d,
                                                      GlobalC::ucell,
                                                      GlobalC::ORB,
                                                      GlobalC::GridD,
                                                      *this->LOWF.ParaV);
                    GlobalC::ld.save_npy_orbital_precalc(GlobalC::ucell.nat,nks);
                    GlobalC::ld.cal_o_delta_k(dm_bandgap_k, *this->LOWF.ParaV, GlobalC::kv.nks);
                    GlobalC::ld.save_npy_o(deepks_bands - GlobalC::ld.o_delta,"o_base.npy",nks);
                } //end deepks_scf multi-k
            } //end deepks_scf == 1
            else //deepks_scf == 0
            { 
                GlobalC::ld.save_npy_o(deepks_bands, "o_base.npy",nks); // no scf, o_tot=o_base
            } //end deepks_scf == 0
        }// end bandgap label
    } // end deepks_out_labels
#endif

    // 3. DeePKS PDM and descriptor
#ifdef __DEEPKS
    const Parallel_Orbitals* pv = this->LOWF.ParaV;
    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        if (GlobalV::GAMMA_ONLY_LOCAL)
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
        GlobalC::ld.check_projected_dm(); //print out the projected dm for NSCF calculaiton
        GlobalC::ld.cal_descriptor(); // final descriptor
        GlobalC::ld.check_descriptor(GlobalC::ucell);

        if (GlobalV::deepks_out_labels)
            GlobalC::ld.save_npy_d(GlobalC::ucell.nat); // libnpy needed
    }

    if (GlobalV::deepks_scf)
    {
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            GlobalC::ld.cal_e_delta_band(this->LOC.dm_gamma, pv->trace_loc_row, pv->trace_loc_col, pv->nrow);
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
        std::cout << "E_delta_band = " << std::setprecision(8) << GlobalC::ld.e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << GlobalC::ld.e_delta_band * ModuleBase::Ry_to_eV << " eV"
                  << std::endl;
        std::cout << "E_delta_NN= " << std::setprecision(8) << GlobalC::ld.E_delta << " Ry"
                  << " = " << std::setprecision(8) << GlobalC::ld.E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
#endif
    // 4. some outputs
#ifdef __EXX
    if(INPUT.rpa)
    {
        // ModuleRPA::DFT_RPA_interface rpa_interface(GlobalC::exx_info.info_global);
        // rpa_interface.rpa_exx_lcao().info.files_abfs = GlobalV::rpa_orbitals;
        // rpa_interface.out_for_RPA(*(this->LOWF.ParaV), *(this->psi), this->LOC, this->pelec);
        RPA_LRI<double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(MPI_COMM_WORLD,this->LOC, *this->LOWF.ParaV);
        rpa_lri_double.init(MPI_COMM_WORLD);
        rpa_lri_double.out_for_RPA(*(this->LOWF.ParaV), *(this->psi), this->LOC, this->pelec);
    }
#endif
    if (hsolver::HSolverLCAO::out_mat_hsR)
    {
        if( !(GlobalV::CALCULATION=="md" && (istep%hsolver::HSolverLCAO::out_hsR_interval!=0)) )
        {
            ModuleIO::output_HS_R(istep, this->pelec->pot->get_effective_v(), this->UHM); // LiuXh add 2019-07-15
        } // LiuXh add 2019-07-15
    }

    if (hsolver::HSolverLCAO::out_mat_t)
    {
        if( !(GlobalV::CALCULATION=="md" && (istep%hsolver::HSolverLCAO::out_hsR_interval!=0)) )
        {
            ModuleIO::output_T_R(istep, this->UHM); // LiuXh add 2019-07-15
        } // LiuXh add 2019-07-15
    }

    if (hsolver::HSolverLCAO::out_mat_dh)
    {
        if( !(GlobalV::CALCULATION=="md" && (istep%hsolver::HSolverLCAO::out_hsR_interval!=0)) )
        {
            ModuleIO::output_dH_R(istep, this->pelec->pot->get_effective_v(), this->UHM); // LiuXh add 2019-07-15
        } // LiuXh add 2019-07-15
    }

    // add by jingan for out r_R matrix 2019.8.14
    if(INPUT.out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(*this->LOWF.ParaV);

        if (hsolver::HSolverLCAO::out_mat_hsR)
        {
            r_matrix.out_rR_other(istep, this->LM.output_R_coor);
        }
        else
        {
            r_matrix.out_rR(istep);
        }
    }

    if(!GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {
        RA.delete_grid();
    }
}

bool ESolver_KS_LCAO::do_after_converge(int& iter)
{
#ifdef __EXX

    // Add EXX operator
    auto add_exx_operator = [&]()
    {
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            hamilt::Operator<double>* exx
                = new hamilt::OperatorEXX<hamilt::OperatorLCAO<double>>(
                    &LM,
                    nullptr, //no explicit call yet
                    &(LM.Hloc)
                );
            p_hamilt->opsd->add(exx);
        }
        else
        {
            hamilt::Operator<std::complex<double>>* exx
                = new hamilt::OperatorEXX<hamilt::OperatorLCAO<std::complex<double>>>(
                    &LM,
                    nullptr, //no explicit call yet
                    &(LM.Hloc2)
                );
            p_hamilt->ops->add(exx);
        }
    };

    if (GlobalC::exx_info.info_global.cal_exx)
    {
        //no separate_loop case
        if (!GlobalC::exx_info.info_global.separate_loop)
        {
            GlobalC::exx_info.info_global.hybrid_step = 1;
                
            //in no_separate_loop case, scf loop only did twice
            //in first scf loop, exx updated once in beginning,
            //in second scf loop, exx updated every iter

            if(this->two_level_step)
            {
                return true;
            }
            else
            {
                //update exx and redo scf
                XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
                iter = 0;
                std::cout << " Entering 2nd SCF, where EXX is updated" << std::endl;
                this->two_level_step++;

                add_exx_operator();

                return false;
            }
        }
        //has separate_loop case
        //exx converged or get max exx steps
        else if(this->two_level_step == GlobalC::exx_info.info_global.hybrid_step || (iter==1 && this->two_level_step!=0))
        {
            return true;
        }
        else
        {
            //update exx and redo scf
            if (two_level_step == 0)
            {
                add_exx_operator();
                XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
            }

            //GlobalC::exx_lcao.cal_exx_elec(this->LOC, this->LOWF.wfc_k_grid);
			if(GlobalC::exx_info.info_ri.real_number)
				GlobalC::exx_lri_double.cal_exx_elec(this->LOC, *this->LOWF.ParaV);
			else
				GlobalC::exx_lri_complex.cal_exx_elec(this->LOC, *this->LOWF.ParaV);
            iter = 0;
            std::cout << " Updating EXX and rerun SCF" << std::endl;
            this->two_level_step++;
            return false;
        }
    }
#endif // __EXX
    return true;
}

} // namespace ModuleESolver
