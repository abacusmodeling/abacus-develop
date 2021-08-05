#include "run_md_pw.h"
#include "forces.h"
#include "stress_pw.h"
#include "global.h" // use chr.
#include "vdwd2.h"
#include "vdwd3.h"
#include "vdwd2_parameters.h"
#include "vdwd3_parameters.h"
#include "pw_complement.h"
#include "pw_basis.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "../module_md/MD_basic.h"

void Run_MD_PW::md_ions_pw(void)
{
    TITLE("Run_MD_PW", "md_ions_pw");
    timer::tick("Run_MD_PW", "md_ions_pw");

    if (GlobalV::OUT_LEVEL == "i")
    {
        std::cout << std::setprecision(12);
        std::cout << " " << std::setw(7) << "ISTEP"
             << std::setw(5) << "NE"
             << std::setw(15) << "ETOT(eV)"
             << std::setw(15) << "EDIFF(eV)"
             << std::setw(15) << "MAX_F(eV/A)"
             << std::setw(15) << "TRADIUS(Bohr)"
             << std::setw(8) << "UPDATE"
             << std::setw(11) << "ETIME(MIN)"
             << std::setw(11) << "FTIME(MIN)"
             << std::endl;
    }

    // allocation for ion movement.
    if (GlobalV::FORCE)
    {
        IMM.allocate();
        CE.allocate_ions();
    }

    if (GlobalV::STRESS) // pengfei Li 2018-05-14
    {
        LCM.allocate();
    }

    MD_basic mdb(INPUT.mdp, GlobalC::ucell);
    int mdtype = INPUT.mdp.mdtype;

    this->istep = 1;
    bool stop = false;

    while (istep <= GlobalV::NSTEP && !stop)
    {
        time_t estart = time(NULL);

        if (GlobalV::OUT_LEVEL == "ie")
        {
            std::cout << " -------------------------------------------" << std::endl;    
            std::cout << " STEP OF MOLECULAR DYNAMICS : " << istep << std::endl;
            std::cout << " -------------------------------------------" << std::endl;
            GlobalV::ofs_running << " -------------------------------------------" << std::endl;
            GlobalV::ofs_running << " STEP OF MOLECULAR DYNAMICS : " << istep << std::endl;
            GlobalV::ofs_running << " -------------------------------------------" << std::endl;
        }

    //----------------------------------------------------------
    // about vdw, jiyy add vdwd3 and linpz add vdwd2
    //----------------------------------------------------------	
        if(INPUT.vdw_method=="d2")
        {
            // setup vdwd2 parameters
	        GlobalC::vdwd2_para.initial_parameters(INPUT);
	        GlobalC::vdwd2_para.initset(GlobalC::ucell);
        }
        if(INPUT.vdw_method=="d3_0" || INPUT.vdw_method=="d3_bj")
        {
            GlobalC::vdwd3_para.initial_parameters(INPUT);
        }
        if (GlobalC::vdwd2_para.flag_vdwd2) //Peize Lin add 2014-04-03, update 2021-03-09
        {
            Vdwd2 vdwd2(GlobalC::ucell, GlobalC::vdwd2_para);
            vdwd2.cal_energy();
            GlobalC::en.evdw = vdwd2.get_energy();
        }
        if (GlobalC::vdwd3_para.flag_vdwd3) //jiyy add 2019-05-18, update 2021-05-02
        {
            Vdwd3 vdwd3(GlobalC::ucell, GlobalC::vdwd3_para);
            vdwd3.cal_energy();
            GlobalC::en.evdw = vdwd3.get_energy();
        }

        // mohan added eiter to count for the electron iteration number, 2021-01-28
        int eiter = 0;
        if (GlobalV::CALCULATION == "md")
        {
#ifdef __LCAO
            if (Exx_Global::Hybrid_Type::No == GlobalC::exx_global.info.hybrid_type)
            {
#endif
                elec.self_consistent(istep - 1);
                eiter = elec.iter;
#ifdef __LCAO
            }
            else if (Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type)
            {
                throw std::invalid_argument(TO_STRING(__FILE__) + TO_STRING(__LINE__));
            }
            else // Peize Lin add 2019-03-09
            {
                if (GlobalC::exx_global.info.separate_loop)
                {
                    for (size_t hybrid_step = 0; hybrid_step != GlobalC::exx_global.info.hybrid_step; ++hybrid_step)
                    {
                        elec.self_consistent(istep - 1);
                        eiter += elec.iter;
                        if (elec.iter == 1 || hybrid_step == GlobalC::exx_global.info.hybrid_step - 1) // exx converge
                            break;
                        GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
                        GlobalC::exx_lip.cal_exx();
                    }
                }
                else
                {
                    elec.self_consistent(istep - 1);
                    eiter += elec.iter;
                    GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
                    elec.self_consistent(istep - 1);
                    eiter += elec.iter;
                }
            }
#endif
        }
        // mohan added 2021-01-28, perform stochastic calculations
        else if (GlobalV::CALCULATION == "md-sto")
        {
            elec_sto.scf_stochastic(istep - 1);
            eiter = elec_sto.iter;
        }

        CE.update_all_pos(GlobalC::ucell);

        if (GlobalC::pot.out_potential == 2)
        {
            std::stringstream ssp;
            std::stringstream ssp_ave;
            ssp << GlobalV::global_out_dir << "ElecStaticPot";
            ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
            GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        time_t eend = time(NULL);
        time_t fstart = time(NULL);

        if (mdtype == 1 || mdtype == 2)
        {
            mdb.runNVT(istep);
        }
        else if (mdtype == 0)
        {
            mdb.runNVE(istep);
        }
        else if (mdtype == -1)
        {
            stop = mdb.runFIRE(istep);
        }
        else
        {
            WARNING_QUIT("opt_ions", "mdtype should be -1~2!");
        }

        time_t fend = time(NULL);

        CE.save_pos_next(GlobalC::ucell);

        //xiaohui add CE.istep = istep 2014-07-07
        CE.update_istep(istep);

        // charge extrapolation if istep>0.
        CE.extrapolate_charge();

        //reset local potential and initial wave function
        GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);
        GlobalV::ofs_running << " Setup the new wave functions?" << std::endl;
        GlobalC::wf.wfcinit();

        if (GlobalV::OUT_LEVEL == "i")
        {
            double etime_min = difftime(eend, estart) / 60.0;
            double ftime_min = difftime(fend, fstart) / 60.0;
            std::stringstream ss;
            ss << GlobalV::MOVE_IONS << istep;

            std::cout << " " << std::setw(7) << ss.str()
                 << std::setw(5) << eiter
                 << std::setw(15) << std::setprecision(6) << GlobalC::en.etot * Ry_to_eV
                 << std::setw(15) << IMM.get_ediff() * Ry_to_eV
                 << std::setprecision(3)
                 << std::setw(15) << IMM.get_largest_grad() * Ry_to_eV / 0.529177
                 << std::setw(15) << IMM.get_trust_radius()
                 << std::setw(8) << IMM.get_update_iter()
                 << std::setprecision(2) << std::setw(11) << etime_min
                 << std::setw(11) << ftime_min << std::endl;
        }

        ++istep;
    }

    if (GlobalV::OUT_LEVEL == "i")
    {
        std::cout << " ION DYNAMICS FINISHED :)" << std::endl;
    }

    timer::tick("Run_MD_PW", "md_ions_pw");
    return;
}

void Run_MD_PW::md_cells_pw()
{
    TITLE("Run_MD_PW", "md_cells_pw");
    timer::tick("Run_MD_PW", "md_cells_pw");

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
    GlobalC::pot.init_pot(0, GlobalC::pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    GlobalC::pot.newd();

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
        throw std::invalid_argument(TO_STRING(__FILE__) + TO_STRING(__LINE__));
    }
#endif

    DONE(GlobalV::ofs_running, "INIT BASIS");

    // ion optimization begins
    // electron density optimization is included in ion optimization

    this->md_ions_pw();

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    timer::tick("Run_MD_PW", "md_cells_pw");
}