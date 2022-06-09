#include "esolver_ks_lcao_tddft.h"

//--------------temporary----------------------------
#include "../module_base/global_function.h"
#include "../src_io/print_info.h"
#include "../src_pw/global.h"
#include "input_update.h"
#include "src_io/chi0_hilbert.h"
#include "src_lcao/ELEC_cbands_gamma.h"
#include "src_lcao/ELEC_cbands_k.h"
#include "src_lcao/ELEC_evolve.h"
#include "src_lcao/dftu.h"
#include "src_lcao/dmft.h"
#include "src_pw/occupy.h"
#include "src_pw/symmetry_rho.h"
#include "src_pw/threshold_elec.h"

#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif
//-----force& stress-------------------
#include "src_lcao/FORCE_STRESS.h"

//---------------------------------------------------

namespace ModuleESolver
{

ESolver_KS_LCAO::ESolver_KS_LCAO_TDDFT()
{
    classname = "ESolver_KS_LCAO_TDDFT";
    basisname = "LCAO";
}
ESolver_KS_LCAO::~ESolver_KS_LCAO_TDDFT()
{
    this->orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
}

void ESolver_KS_LCAO_TDDFT::Init(Input& inp, UnitCell_pseudo& ucell)
{
    // setup GlobalV::NBANDS
    // Yu Liu add 2021-07-03
    GlobalC::CHR.cal_nelec();

    // it has been established that that
    // xc_func is same for all elements, therefore
    // only the first one if used
    if (ucell.atoms[0].xc_func == "HSE" || ucell.atoms[0].xc_func == "PBE0")
    {
        XC_Functional::set_xc_type("pbe");
    }
    else
    {
        XC_Functional::set_xc_type(ucell.atoms[0].xc_func);
    }

    // ucell.setup_cell( GlobalV::global_pseudo_dir , GlobalV::stru_file , GlobalV::ofs_running, GlobalV::NLOCAL,
    // GlobalV::NBANDS);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set(GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(ucell, GlobalC::kv);

    //--------------------------------------
    // cell relaxation should begin here
    //--------------------------------------

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx << " * " << GlobalC::pw.ny << " * " << GlobalC::pw.nz
              << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx << " * " << GlobalC::pw.nby << " * " << GlobalC::pw.nbz
              << std::endl;

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx,
                        GlobalC::pw.ncy,
                        GlobalC::pw.ncz,
                        GlobalC::pw.nczp,
                        GlobalC::pw.nrxx,
                        GlobalC::pw.nbz,
                        GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04
    // Calculate Structure factor
    GlobalC::pw.setup_structure_factor();

    // Inititlize the charge density.
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT CHARGE");

    // Initializee the potential.
    GlobalC::pot.allocate(GlobalC::pw.nrxx);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

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
    int ion_step = 0;
    GlobalC::pot.init_pot(ion_step, GlobalC::pw.strucFac);

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

    // init Psi
    if (GlobalV::GAMMA_ONLY_LOCAL)
        this->LOWF.wfc_gamma.resize(GlobalV::NSPIN);
    else
    {
        this->LOWF.wfc_k.resize(GlobalC::kv.nks);
        this->LOWF.wfc_k_laststep.resize(GlobalC::kv.nks);
    }
    // pass Hamilt-pointer to Operator
    this->UHM.genH.LM = this->UHM.LM = &this->LM;
    // pass basis-pointer to EState and Psi
    this->LOC.ParaV = this->LOWF.ParaV = this->LM.ParaV;
}

void ESolver_KS_LCAO_TDDFT::eachiterinit(const int istep, const int iter)
{
    std::string ufile = "CHANGE";
    Update_input UI;
    UI.init(ufile);

    if (INPUT.dft_plus_u)
        GlobalC::dftu.iter_dftu = iter;

    // mohan add 2010-07-16
    // used for pulay mixing.
    if (iter == 1)
    {
        GlobalC::CHR.set_new_e_iteration(true);
    }
    else
    {
        GlobalC::CHR.set_new_e_iteration(false);
    }

    if (GlobalV::FINAL_SCF && iter == 1)
    {
        GlobalC::CHR.irstep = 0;
        GlobalC::CHR.idstep = 0;
        GlobalC::CHR.totstep = 0;
    }

    // mohan update 2012-06-05
    GlobalC::en.calculate_harris(1);

    // mohan move it outside 2011-01-13
    // first need to calculate the weight according to
    // electrons number.
    // mohan add iter > 1 on 2011-04-02
    // because the GlobalC::en.ekb has not value now.
    // so the smearing can not be done.
    if (iter > 1 && istep <= 1 && GlobalV::ocp == 0)
        Occupy::calculate_weights();

    if (istep <= 1 && GlobalV::ocp == 1)
    {
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                GlobalC::wf.wg(ik, ib) = GlobalV::ocp_kb[ik * GlobalV::NBANDS + ib];
            }
        }
    }

    if (GlobalC::wf.init_wfc == "file")
    {
        if (iter == 1)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // The occupation should be read in together.
            // Occupy::calculate_weights(); //mohan add 2012-02-15

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.
            this->LOC.sum_bands(this->UHM);
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
            GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
            GlobalC::en.delta_escf();
            if (ELEC_evolve::td_vext == 0)
            {
                GlobalC::pot.set_vr_eff();
            }
            else
            {
                GlobalC::pot.set_vrs_tddft(istep);
            }
        }
    }
}

void ESolver_KS_LCAO_TDDFT::hamilt2density(int istep, int iter, double ethr)
{
    // (1) calculate the bands.
    // mohan add 2021-02-09
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        ELEC_cbands_gamma::cal_bands(istep, this->UHM, this->LOWF, this->LOC.dm_gamma);
    }
    else
    {
        if (istep >= 2)
        {
            ELEC_evolve::evolve_psi(istep, this->UHM, this->LOWF);
        }
        else
        {
            ELEC_cbands_k::cal_bands(istep, this->UHM, this->LOWF, this->LOC.dm_k);
        }
    }
    //-----------------------------------------------------------
    // only deal with charge density after both wavefunctions.
    // are calculated.
    //-----------------------------------------------------------
    if (GlobalV::GAMMA_ONLY_LOCAL && GlobalV::NSPIN == 2 && GlobalV::CURRENT_SPIN == 0)
        return;

    GlobalC::en.eband = 0.0;
    GlobalC::en.ef = 0.0;
    GlobalC::en.ef_up = 0.0;
    GlobalC::en.ef_dw = 0.0;

    // demet is included into eband.
    // if(GlobalV::DIAGO_TYPE!="selinv")
    {
        GlobalC::en.demet = 0.0;
    }

    // (2)
    GlobalC::CHR.save_rho_before_sum_band();

    // (3) sum bands to calculate charge density
    if (istep <= 1 && istep <= 1 && GlobalV::ocp == 0)
        Occupy::calculate_weights();

    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        GlobalC::en.print_band(ik);
    }

    // if selinv is used, we need this to calculate the charge
    // using density matrix.
    this->LOC.sum_bands(this->UHM);

    // (4) mohan add 2010-06-24
    // using new charge density.
    GlobalC::en.calculate_harris(2);

    // (5) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, GlobalC::CHR, GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
    }

    // (6) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization();

    // resume codes!
    //-------------------------------------------------------------------------
    // this->GlobalC::LOWF.init_Cij( 0 ); // check the orthogonality of local orbital.
    // GlobalC::CHR.sum_band(); use local orbital in plane wave basis to calculate bands.
    // but must has evc first!
    //-------------------------------------------------------------------------

    // (7) calculate delta energy
    GlobalC::en.deband = GlobalC::en.delta_e();

    // (8) store wfc
#ifdef __MPI
    const Parallel_Orbitals* pv = this->LOWF.ParaV;
    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        this->LOWF.wfc_k_laststep[ik].create(pv->ncol_bands, pv->nrow);
        this->LOWF.wfc_k_laststep[ik] = this->LOWF.wfc_k[ik];
    }

    if (istep > 1)
        this->cal_edm_tddft();
#else

    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        this->LOWF.wfc_k_laststep[ik].create(GlobalV::NBANDS, GlobalV::NLOCAL);
        this->LOWF.wfc_k_laststep[ik] = this->LOWF.wfc_k[ik];
    }
    /*
        this->LM_md.Hloc2_laststep.resize(GlobalV::NLOCAL*GlobalV::NLOCAL);
    ModuleBase::GlobalFunc::ZEROS(this->LM_md.Hloc2_laststep.data(),GlobalV::NLOCAL*GlobalV::NLOCAL);
        for (int i =0 ; i< GlobalV::NLOCAL*GlobalV::NLOCAL;++i)
           {
               this->LM_md.Hloc2_laststep[i]=this->LM_md.Hloc2[i];
           }
    */
    if (istep > 1)
        this->cal_edm_tddft();
#endif
}

void ESolver_KS_LCAO_TDDFT::updatepot(const int istep, const int iter)
{
    // (9) Calculate new potential according to new Charge Density.

    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (GlobalC::pot.out_pot < 0) // mohan add 2011-10-10
        {
            GlobalC::pot.out_pot = -2;
        }
    }
    if (!this->conv_elec)
    {
        GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        GlobalC::en.delta_escf();
    }
    else
    {
        GlobalC::pot.vnew = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        //(used later for scf correction to the forces )
        GlobalC::pot.vnew -= GlobalC::pot.vr;
        GlobalC::en.descf = 0.0;
    }

    // add Vloc to Vhxc.
    if (ELEC_evolve::td_vext == 0)
    {
        GlobalC::pot.set_vr_eff();
    }
    else
    {
        GlobalC::pot.set_vrs_tddft(istep);
    }
}

} // namespace ModuleESolver

// use the original formula (Hamiltonian matrix) to calculate energy density matrix
void ESolver_KS_LCAO_TDDFT::cal_edm_tddft()
{
    this->LOC.edm_k_tddft.resize(GlobalC::kv.nks);
    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
#ifdef __MPI
        this->LOC.edm_k_tddft[ik].create(this->LOC.ParaV->ncol, this->LOC.ParaV->nrow);
        complex<double>* Htmp = new complex<double>[this->LM->ParaV->nloc];
        complex<double>* Sinv = new complex<double>[this->LM->ParaV->nloc];
        complex<double>* tmp1 = new complex<double>[this->LM->ParaV->nloc];
        complex<double>* tmp2 = new complex<double>[this->LM->ParaV->nloc];
        complex<double>* tmp3 = new complex<double>[this->LM->ParaV->nloc];
        complex<double>* tmp4 = new complex<double>[this->LM->ParaV->nloc];
        complex<double>* tmp5 = new complex<double>[this->LM->ParaV->nloc];
        ModuleBase::GlobalFunc::ZEROS(Htmp, this->LM->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(Sinv, this->LM->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp1, this->LM->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp2, this->LM->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp3, this->LM->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp4, this->LM->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp5, this->LM->ParaV->nloc);
        const int inc = 1;
        int nrow = this->LM->ParaV->nrow;
        int ncol = this->LM->ParaV->ncol;
        zcopy_(&this->LM->ParaV->nloc, this->LM->Hloc2.data(), &inc, Htmp, &inc);
        zcopy_(&this->LM->ParaV->nloc, this->LM->Sloc2.data(), &inc, Sinv, &inc);

        int* ipiv = new int[this->LM->ParaV->nloc];
        int info;
        const int one_int = 1;
        pzgetrf_(&GlobalV::NLOCAL, &GlobalV::NLOCAL, Sinv, &one_int, &one_int, this->LM->ParaV->desc, ipiv, &info);

        int LWORK = -1, liWORK = -1;
        std::vector<std::complex<double>> WORK(1, 0);
        std::vector<int> iWORK(1, 0);

        pzgetri_(&GlobalV::NLOCAL,
                 Sinv,
                 &one_int,
                 &one_int,
                 this->LM->ParaV->desc,
                 ipiv,
                 WORK.data(),
                 &LWORK,
                 iWORK.data(),
                 &liWORK,
                 &info);

        LWORK = WORK[0].real();
        WORK.resize(LWORK, 0);
        liWORK = iWORK[0];
        iWORK.resize(liWORK, 0);

        pzgetri_(&GlobalV::NLOCAL,
                 Sinv,
                 &one_int,
                 &one_int,
                 this->LM->ParaV->desc,
                 ipiv,
                 WORK.data(),
                 &LWORK,
                 iWORK.data(),
                 &liWORK,
                 &info);

        const char N_char = 'N', T_char = 'T';
        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};
        const complex<double> half_float[2] = {0.5, 0.0};
        pzgemm_(&T_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                this->LOC.dm_k[ik].c,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                Htmp,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                &zero_float[0],
                tmp1,
                &one_int,
                &one_int,
                this->LM->ParaV->desc);

        pzgemm_(&N_char,
                &N_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                tmp1,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                Sinv,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                &zero_float[0],
                tmp2,
                &one_int,
                &one_int,
                this->LM->ParaV->desc);

        pzgemm_(&N_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                Sinv,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                Htmp,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                &zero_float[0],
                tmp3,
                &one_int,
                &one_int,
                this->LM->ParaV->desc);

        pzgemm_(&N_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                tmp3,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                this->LOC.dm_k[ik].c,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                &zero_float[0],
                tmp4,
                &one_int,
                &one_int,
                this->LM->ParaV->desc);

        pzgeadd_(&N_char,
                 &GlobalV::NLOCAL,
                 &GlobalV::NLOCAL,
                 &half_float[0],
                 tmp2,
                 &one_int,
                 &one_int,
                 this->LM->ParaV->desc,
                 &half_float[0],
                 tmp4,
                 &one_int,
                 &one_int,
                 this->LM->ParaV->desc);

        pztranu_(&GlobalV::NLOCAL,
                 &GlobalV::NLOCAL,
                 &one_float[0],
                 tmp4,
                 &one_int,
                 &one_int,
                 this->LM->ParaV->desc,
                 &zero_float[0],
                 tmp5,
                 &one_int,
                 &one_int,
                 this->LM->ParaV->desc);
        zcopy_(&this->LM->ParaV->nloc, tmp5, &inc, this->LOC.edm_k_tddft[ik].c, &inc);

        delete[] Htmp;
        delete[] Sinv;
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] tmp4;
        delete[] tmp5;
#else
        this->LOC.edm_k_tddft[ik].create(this->LOC.ParaV->ncol, this->LOC.ParaV->nrow);
        ModuleBase::ComplexMatrix Sinv(GlobalV::NLOCAL, GlobalV::NLOCAL);
        ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                Htmp(i, j) = this->LM->Hloc2[i * GlobalV::NLOCAL + j];
                Sinv(i, j) = this->LM->Sloc2[i * GlobalV::NLOCAL + j];
            }
        }
        int INFO;

        int LWORK = 3 * GlobalV::NLOCAL - 1; // tmp
        std::complex<double>* WORK = new std::complex<double>[LWORK];
        ModuleBase::GlobalFunc::ZEROS(WORK, LWORK);
        int IPIV[GlobalV::NLOCAL];

        LapackConnector::zgetrf(GlobalV::NLOCAL, GlobalV::NLOCAL, Sinv, GlobalV::NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri(GlobalV::NLOCAL, Sinv, GlobalV::NLOCAL, IPIV, WORK, LWORK, &INFO);

        this->LOC.edm_k_tddft[ik] = 0.5 * (Sinv * Htmp * this->LOC.dm_k[ik] + this->LOC.dm_k[ik] * Htmp * Sinv);
#endif
    }
    return;
}