#include "esolver_ks.h"
#include <iostream>
#include "time.h"
#include "../module_io/print_info.h"
#ifdef __MPI
#include "mpi.h"
#else
#include "chrono"
#endif

//--------------Temporary----------------
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_base/timer.h"
//---------------------------------------

namespace ModuleESolver
{

    template<typename FPTYPE, typename Device>
    ESolver_KS<FPTYPE, Device>::ESolver_KS()
    {
        classname = "ESolver_KS";
        basisname = "PLEASE ADD BASISNAME FOR CURRENT ESOLVER.";
        scf_thr = GlobalV::SCF_THR;
        drho = 0.0;
        maxniter = GlobalV::SCF_NMAX;
        niter = maxniter;
        out_freq_elec = GlobalV::OUT_FREQ_ELEC;

        // pw_rho = new ModuleBase::PW_Basis();
        //temporary, it will be removed
        pw_wfc = new ModulePW::PW_Basis_K_Big(GlobalV::device_flag, GlobalV::precision_flag);
        GlobalC::wfcpw = this->pw_wfc; //Temporary
        ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(pw_wfc);
        tmp->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    }

    template<typename FPTYPE, typename Device>
    ESolver_KS<FPTYPE, Device>::~ESolver_KS()
    {
        delete this->pw_wfc;
        delete this->p_hamilt;
        delete this->phsol;
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::Init(Input& inp, UnitCell& ucell)
    {
        ESolver_FP::Init(inp,ucell);
        chr.cal_nelec();

        /* it has been established that that
         xc_func is same for all elements, therefore
         only the first one if used*/
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

        // symmetry analysis should be performed every time the cell is changed
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
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

        if(GlobalV::BASIS_TYPE=="pw" || GlobalV::CALCULATION=="ienvelope")
        {
            //Envelope function is calculated as lcao_in_pw
            //new plane wave basis
    #ifdef __MPI
            this->pw_wfc->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
    #endif
            this->pw_wfc->initgrids(ucell.lat0, ucell.latvec, GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz);
            this->pw_wfc->initparameters(false, inp.ecutwfc, GlobalC::kv.nks, GlobalC::kv.kvec_d.data());
    #ifdef __MPI
            if(INPUT.pw_seed > 0)    MPI_Allreduce(MPI_IN_PLACE, &this->pw_wfc->ggecut, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
            //qianrui add 2021-8-13 to make different kpar parameters can get the same results
    #endif
            this->pw_wfc->setuptransform();
            for(int ik = 0 ; ik < GlobalC::kv.nks; ++ik)   GlobalC::kv.ngk[ik] = this->pw_wfc->npwk[ik];
            this->pw_wfc->collect_local_pw(); 
            this->print_wfcfft(inp, GlobalV::ofs_running);
        }
        // initialize the real-space uniform grid for FFT and parallel
        // distribution of plane waves
        GlobalC::Pgrid.init(GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz, GlobalC::rhopw->nplane,
            GlobalC::rhopw->nrxx, GlobalC::bigpw->nbz, GlobalC::bigpw->bz); // mohan add 2010-07-22, update 2011-05-04
            
        // Calculate Structure factor
        GlobalC::sf.setup_structure_factor(&GlobalC::ucell, GlobalC::rhopw);

        // Initialize charge extrapolation
        CE.Init_CE();
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::hamilt2density(const int istep, const int iter, const FPTYPE ethr)
    {
        ModuleBase::timer::tick(this->classname, "hamilt2density");
        //Temporarily, before HSolver is constructed, it should be overrided by
        //LCAO, PW, SDFT and TDDFT.
        //After HSolver is constructed, LCAO, PW, SDFT should delete their own
        //hamilt2density() and use:
        //this->phsol->solve(this->phamilt, this->pes, this->wf, ETHR);
        ModuleBase::timer::tick(this->classname, "hamilt2density");
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::print_wfcfft(Input& inp, ofstream &ofs)
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
        FPTYPE ecut = INPUT.ecutwfc;
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

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::Run(const int istep, UnitCell& ucell)
    {
        if (!(GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
            || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax"))
        {
            this->othercalculation(istep);
        }
        else
        {
            ModuleBase::timer::tick(this->classname, "Run");

            this->beforescf(istep); //Something else to do before the iter loop
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");
            if(this->maxniter > 0)  this->printhead(); //print the headline on the screen.

            bool firstscf = true;
            this->conv_elec = false;
            this->niter = this->maxniter;
            for (int iter = 1; iter <= this->maxniter; ++iter)
            {
                writehead(GlobalV::ofs_running, istep, iter);
#ifdef __MPI
                auto iterstart = MPI_Wtime();
#else
                auto iterstart = std::chrono::system_clock::now();
#endif
                FPTYPE diag_ethr = this->phsol->set_diagethr(istep, iter, drho);
                eachiterinit(istep, iter);
                this->hamilt2density(istep, iter, diag_ethr);
                
                //<Temporary> It may be changed when more clever parallel algorithm is put forward.
                //When parallel algorithm for bands are adopted. Density will only be treated in the first group.
                //(Different ranks should have abtained the same, but small differences always exist in practice.)
                //Maybe in the future, density and wavefunctions should use different parallel algorithms, in which 
                //they do not occupy all processors, for example wavefunctions uses 20 processors while density uses 10.
                if(GlobalV::MY_STOGROUP == 0)
                {
                    // FPTYPE drho = this->estate.caldr2(); 
                    // EState should be used after it is constructed.

                    drho = GlobalC::CHR_MIX.get_drho(pelec->charge, GlobalV::nelec);
                    FPTYPE hsolver_error = 0.0;
                    if (firstscf)
                    {
                        firstscf = false;
                        hsolver_error = this->phsol->cal_hsolerror();
                        // The error of HSolver is larger than drho, so a more precise HSolver should be excuconv_elected.
                        if (hsolver_error > drho)
                        {
                            diag_ethr = this->phsol->reset_diagethr(GlobalV::ofs_running, hsolver_error, drho);
                            this->hamilt2density(istep, iter, diag_ethr);
                            drho = GlobalC::CHR_MIX.get_drho(pelec->charge, GlobalV::nelec);
                            hsolver_error = this->phsol->cal_hsolerror();
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
                        GlobalC::CHR_MIX.mix_rho(iter, pelec->charge);
                    }
                }
#ifdef __MPI
		        MPI_Bcast(&drho, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		        MPI_Bcast(&this->conv_elec, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		        MPI_Bcast(pelec->charge->rho[0], GlobalC::rhopw->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

                // Hamilt should be used after it is constructed.
                // this->phamilt->update(conv_elec);
                updatepot(istep, iter);
                eachiterfinish(iter);
#ifdef __MPI
                FPTYPE duration = (FPTYPE)(MPI_Wtime() - iterstart);
#else
                FPTYPE duration = (std::chrono::system_clock::now() - iterstart).count() / CLOCKS_PER_SEC;
#endif
                printiter(iter, drho, duration, diag_ethr);
                if (this->conv_elec)
                {
                    this->niter = iter;
                    bool stop = this->do_after_converge(iter);
                    if(stop) break;
                }
            }
            afterscf(istep);

            ModuleBase::timer::tick(this->classname, "Run");
        }       

        return;
    };

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::printhead()
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

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::printiter(const int iter, const FPTYPE drho, const FPTYPE duration, const FPTYPE ethr)
    {
        GlobalC::en.print_etot(this->conv_elec, iter, drho, duration, ethr);
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS<FPTYPE, Device>::writehead(std::ofstream& ofs_running, const int istep, const int iter)
    {
        ofs_running
            << "\n "
            << this->basisname
            << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
            << "  ELEC=" << std::setw(4) << iter
            << "--------------------------------\n";
    }

    template<typename FPTYPE, typename Device>
    int ESolver_KS<FPTYPE, Device>::getniter()
    {
        return this->niter;
    }

template class ESolver_KS<float, psi::DEVICE_CPU>;
template class ESolver_KS<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS<float, psi::DEVICE_GPU>;
template class ESolver_KS<double, psi::DEVICE_GPU>;
#endif
}