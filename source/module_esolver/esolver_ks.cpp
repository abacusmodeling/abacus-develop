#include "esolver_ks.h"
#include <iostream>
#include <algorithm>
#include "time.h"
#include "../src_io/print_info.h"
#ifdef __MPI
#include "mpi.h"
#endif

//--------------Temporary----------------
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_pw/charge_broyden.h"
#include "../module_base/timer.h"
//---------------------------------------

namespace ModuleESolver
{

    ESolver_KS::ESolver_KS()
    {
        classname = "ESolver_KS";
        basisname = "PLEASE ADD BASISNAME FOR CURRENT ESOLVER.";
        diag_ethr = GlobalV::PW_DIAG_THR;
        scf_thr = GlobalV::SCF_THR;
        drho = 0.0;
        maxniter = GlobalV::SCF_NMAX;
        niter = maxniter;
        out_freq_elec = GlobalV::OUT_FREQ_ELEC;

        // pw_rho = new ModuleBase::PW_Basis();
        //temporary, it will be removed
        pw_wfc = new ModulePW::PW_Basis_K_Big(); 
        GlobalC::wfcpw = this->pw_wfc; //Temporary
        ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(pw_wfc);
        tmp->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    }

    void ESolver_KS::Init(Input& inp, UnitCell_pseudo& ucell)
    {
        ESolver_FP::Init(inp,ucell);
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

        //new plane wave basis
        this->pw_wfc->initgrids(ucell.lat0, ucell.latvec, GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz,
                                GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
        this->pw_wfc->initparameters(false, inp.ecutwfc, GlobalC::kv.nks, GlobalC::kv.kvec_d.data());
#ifdef __MPI
        if(INPUT.pw_seed > 0)    MPI_Allreduce(MPI_IN_PLACE, &this->pw_wfc->ggecut, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
        //qianrui add 2021-8-13 to make different kpar parameters can get the same results
#endif
        this->pw_wfc->setuptransform();
        for(int ik = 0 ; ik < GlobalC::kv.nks; ++ik)   GlobalC::kv.ngk[ik] = this->pw_wfc->npwk[ik];
        this->pw_wfc->collect_local_pw(); 
        print_wfcfft(inp, GlobalV::ofs_running);

        // initialize the real-space uniform grid for FFT and parallel
        // distribution of plane waves
        GlobalC::Pgrid.init(GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz, GlobalC::rhopw->nplane,
            GlobalC::rhopw->nrxx, GlobalC::bigpw->nbz, GlobalC::bigpw->bz); // mohan add 2010-07-22, update 2011-05-04
        // Calculate Structure factor
        GlobalC::sf.setup_structure_factor(&GlobalC::ucell, GlobalC::rhopw);

        // Inititlize the charge density.
        GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT CHARGE");
        // Initializee the potential.
        GlobalC::pot.allocate(GlobalC::rhopw->nrxx);
    }

    void ESolver_KS::hamilt2density(const int istep, const int iter, const double ethr)
    {
        ModuleBase::timer::tick(this->classname, "hamilt2density");
        //Temporarily, before HSolver is constructed, it should be overrided by
        //LCAO, PW, SDFT and TDDFT.
        //After HSolver is constructed, LCAO, PW, SDFT should delete their own
        //hamilt2density() and use:
        //this->phsol->solve(this->phamilt, this->pes, this->wf, ETHR);
        ModuleBase::timer::tick(this->classname, "hamilt2density");
    }

    void ESolver_KS::print_wfcfft(Input& inp, ofstream &ofs)
    {
        ofs << "\n\n\n\n";
	    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	    ofs << " |                                                                    |" << std::endl;
	    ofs << " | Setup plane waves of wave functions:                               |" << std::endl;
	    ofs << " | Use the energy cutoff and the lattice vectors to generate the      |" << std::endl;
	    ofs << " | dimensions of FFT grid. The number of FFT grid on each processor   |" << std::endl;
	    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space is   |" << std::endl;
	    ofs << " | different for charege/potential and wave functions. We also set    |" << std::endl;
	    ofs << " | the 'sticks' for the parallel of FFT. The number of plane wave of  |" << std::endl;
	    ofs << " | each k-point is 'npwk[ik]' in each processor                       |" << std::endl;
	    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	    ofs << "\n\n\n\n";
        ofs << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;
        double ecut = INPUT.ecutwfc;
        if(abs(ecut-this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2) > 1e-6)
        {
            ecut = this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2;
            ofs<<"Energy cutoff for wavefunc is incompatible with nx, ny, nz and it will be reduced!"<<std::endl;

        }  
        ModuleBase::GlobalFunc::OUT(ofs,"energy cutoff for wavefunc (unit:Ry)", ecut);
        ModuleBase::GlobalFunc::OUT(ofs,"fft grid for wave functions", this->pw_wfc->nx,this->pw_wfc->ny,this->pw_wfc->nz);
        ModuleBase::GlobalFunc::OUT(ofs,"number of plane waves",this->pw_wfc->npwtot);
	    ModuleBase::GlobalFunc::OUT(ofs,"number of sticks", this->pw_wfc->nstot);

        ofs << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
        ofs <<" "<< std::setw(8)  << "PROC"<< std::setw(15) << "COLUMNS(POT)"<< std::setw(15) << "PW" << std::endl;
        for (int i = 0; i < GlobalV::NPROC_IN_POOL ; ++i)
        {
            ofs <<" "<<std::setw(8)<< i+1 << std::setw(15) << this->pw_wfc->nst_per[i] << std::setw(15) << this->pw_wfc->npw_per[i] << std::endl;
        }
        ofs << " --------------- sum -------------------" << std::endl;
        ofs << " " << std::setw(8)  << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_wfc->nstot << std::setw(15) << this->pw_wfc->npwtot << std::endl;
        ModuleBase::GlobalFunc::DONE(ofs, "INIT PLANEWAVE");
    }


    void ESolver_KS::Run(const int istep, UnitCell_pseudo& ucell)
    {
        if (!(GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
            || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax" || GlobalV::CALCULATION.substr(0,3) == "sto")
#ifdef __MPI
            || Exx_Global::Hybrid_Type::No != GlobalC::exx_global.info.hybrid_type
#endif
            )
        {
            this->othercalculation(istep);
        }
        else
        {
            ModuleBase::timer::tick(this->classname, "Run");

            this->printhead(); //print the headline on the screen.
            this->beforescf(istep); //Something else to do before the iter loop

            bool firstscf = true;
            this->conv_elec = false;
            this->niter = this->maxniter;
            for (int iter = 1; iter <= this->maxniter; ++iter)
            {
                writehead(GlobalV::ofs_running, istep, iter);
                clock_t iterstart, iterend;
                iterstart = std::clock();
                set_ethr(istep, iter);
                eachiterinit(istep, iter);
                this->hamilt2density(istep, iter, this->diag_ethr);
                
                //<Temporary> It may be changed when more clever parallel algorithm is put forward.
                //When parallel algorithm for bands are adopted. Density will only be treated in the first group.
                //(Different ranks should have abtained the same, but small differences always exist in practice.)
                //Maybe in the future, density and wavefunctions should use different parallel algorithms, in which 
                //they do not occupy all processors, for example wavefunctions uses 20 processors while density uses 10.
                if(GlobalV::MY_STOGROUP == 0)
                {
                // double drho = this->estate.caldr2(); 
                // EState should be used after it is constructed.
                drho = GlobalC::CHR.get_drho();
                double hsolver_error = 0.0;
                if (firstscf)
                {
                    firstscf = false;
                    hsolver_error = this->diag_ethr * std::max(1.0, GlobalC::CHR.nelec);
                    // The error of HSolver is larger than drho, so a more precise HSolver should be excuconv_elected.
                    if (hsolver_error > drho)
                    {
                        reset_diagethr(GlobalV::ofs_running, hsolver_error);
                        this->hamilt2density(istep, iter, this->diag_ethr);
                        drho = GlobalC::CHR.get_drho();
                        hsolver_error = this->diag_ethr * std::max(1.0, GlobalC::CHR.nelec);
                    }
                }

                this->conv_elec = (drho < this->scf_thr);

                // If drho < hsolver_error in the first iter or drho < scf_thr, we do not change rho.
                if (drho < hsolver_error || this->conv_elec)
                {
                    if (drho < hsolver_error)    GlobalV::ofs_warning << " drho < hsolver_error, keep charge density unchanged." << std::endl;
                }
                else
                {
                    //charge mixing
                    //conv_elec = this->estate.mix_rho();
                    GlobalC::CHR.mix_rho(iter);
                }
                }
#ifdef __MPI
		        MPI_Bcast(&drho, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		        MPI_Bcast(&this->conv_elec, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		        MPI_Bcast(GlobalC::CHR.rho[0], GlobalC::rhopw->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

                // Hamilt should be used after it is constructed.
                // this->phamilt->update(conv_elec);
                updatepot(istep, iter);
                eachiterfinish(iter);
                iterend = std::clock();
                double duration = double(iterend - iterstart) / CLOCKS_PER_SEC;
                printiter(iter, drho, duration, diag_ethr);
                if (this->conv_elec)
                {
                    this->niter = iter;
                    break;
                }
            }
            afterscf();

            ModuleBase::timer::tick(this->classname, "Run");
        }       

        return;
    };

    //<Temporary> It should be a function of Diag_H class in the future.
    void ESolver_KS::set_ethr(const int istep, const int iter)
    {
        //It is too complex now and should be modified.
        if (iter == 1)
        {
            if (abs(this->diag_ethr - 1.0e-2) < 1.0e-10)
            {
                if (GlobalC::pot.init_chg == "file")
                {
                    //======================================================
                    // if you think that the starting potential is good
                    // do not spoil it with a louly first diagonalization:
                    // set a strict this->diag_ethr in the input file ()diago_the_init
                    //======================================================
                    this->diag_ethr = 1.0e-5;
                }
                else
                {
                    //=======================================================
                    // starting atomic potential is probably far from scf
                    // don't waste iterations in the first diagonalization
                    //=======================================================
                    this->diag_ethr = 1.0e-2;
                }
            }
            if (GlobalV::FINAL_SCF) this->diag_ethr = 1.0e-2;

            if (GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
            {
                this->diag_ethr = std::max(this->diag_ethr, INPUT.pw_diag_thr);
            }

        }
        else
        {
            if (iter == 2)
            {
                this->diag_ethr = 1.e-2;
            }
            this->diag_ethr = std::min(this->diag_ethr, 0.1 * this->drho / std::max(1.0, GlobalC::CHR.nelec));

        }
        if (GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw"|| GlobalV::CALCULATION.substr(0,3)=="sto")
        {
            this->diag_ethr = 0.0;
        }
    }

    void ESolver_KS::printhead()
    {
        std::cout << " " << std::setw(7) << "ITER";
        if (GlobalV::NSPIN == 2)
        {
            std::cout << std::setw(10) << "TMAG";
            std::cout << std::setw(10) << "AMAG";
        }
        std::cout << std::setw(15) << "ETOT(eV)";
        std::cout << std::setw(15) << "EDIFF(eV)";
        std::cout << std::setw(11) << "DRHO";
        std::cout << std::setw(11) << "TIME(s)" << std::endl;
    }

    void ESolver_KS::printiter(const int iter, const double drho, const double duration, const double ethr)
    {
        GlobalC::en.print_etot(this->conv_elec, iter, drho, duration, ethr);
    }

    void ESolver_KS::writehead(std::ofstream& ofs_running, const int istep, const int iter)
    {
        ofs_running
            << "\n "
            << this->basisname
            << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
            << "  ELEC=" << std::setw(4) << iter
            << "--------------------------------\n";
    }

    void ESolver_KS::reset_diagethr(std::ofstream& ofs_running, const double hsover_error)
    {
        ofs_running << " Notice: Threshold on eigenvalues was too large.\n";
        ModuleBase::WARNING("scf", "Threshold on eigenvalues was too large.");
        ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho << std::endl;
        ofs_running << " Origin diag_ethr = " << this->diag_ethr << std::endl;
        this->diag_ethr = 0.1 * drho / GlobalC::CHR.nelec;
        ofs_running << " New    diag_ethr = " << this->diag_ethr << std::endl;
    }

    int ESolver_KS::getniter()
    {
        return this->niter;
    }


}
