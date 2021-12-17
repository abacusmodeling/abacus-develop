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
#include "../module_md/MD_func.h"
#include "../module_md/FIRE.h"
#include "../module_md/NVE.h"
#include "../module_md/MSST.h"
#include "../module_md/NVT_ADS.h"
#include "../module_md/NVT_NHC.h"
#include "../module_md/Langevin.h"
#include "../src_io/print_info.h"

Run_MD_PW::Run_MD_PW()
{
    cellchange = false;
}

Run_MD_PW::~Run_MD_PW(){}

void Run_MD_PW::md_ions_pw(void)
{
    ModuleBase::TITLE("Run_MD_PW", "md_ions_pw");
    ModuleBase::timer::tick("Run_MD_PW", "md_ions_pw");

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
    CE.allocate_ions();

    // determine the mdtype
    Verlet *verlet;
    if(INPUT.mdp.mdtype == -1)
    {
        verlet = new FIRE(INPUT.mdp, GlobalC::ucell); 
    }
    else if(INPUT.mdp.mdtype == 0)
    {
        verlet = new NVE(INPUT.mdp, GlobalC::ucell); 
    }
    else if(INPUT.mdp.mdtype==1)
    {
        verlet = new NVT_ADS(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.mdtype==2)
    {
        verlet = new NVT_NHC(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.mdtype == 3)
    {
        verlet = new Langevin(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.mdtype==4)
    {
        verlet = new MSST(INPUT.mdp, GlobalC::ucell); 
        cellchange = true;
    }

    // md cycle
    while ( (verlet->step_ + verlet->step_rst_) <= GlobalV::NSTEP && !verlet->stop)
    {
        if(verlet->step_ == 0)
        {
            verlet->setup();
        }
        else
        {
            CE.update_all_pos(GlobalC::ucell);

            verlet->first_half();

            if(cellchange)
            {
                CE.update_istep(1);
            }
            else
            {
                CE.update_istep(verlet->step_);
            }

            CE.extrapolate_charge();
            CE.save_pos_next(GlobalC::ucell);

            if(cellchange)
            {
                Variable_Cell::init_after_vc();
            }

            // reset local potential and initial wave function
            GlobalC::pot.init_pot(verlet->step_, GlobalC::pw.strucFac);
            
            // new wave functions
            GlobalC::wf.wfcinit();

            // update force and virial due to the update of atom positions
            MD_func::force_virial(verlet->step_, verlet->mdp, verlet->ucell, verlet->potential, verlet->force, verlet->virial);

            verlet->second_half();

            MD_func::kinetic_stress(verlet->ucell, verlet->vel, verlet->allmass, verlet->kinetic, verlet->stress);

            verlet->stress += verlet->virial;
        }

        if((verlet->step_ + verlet->step_rst_) % verlet->mdp.dumpfreq == 0)
        {
            Print_Info::print_screen(0, 0, verlet->step_ + verlet->step_rst_);
            verlet->outputMD();

            MD_func::MDdump(verlet->step_ + verlet->step_rst_, verlet->ucell, verlet->virial, verlet->force);
        }

        if((verlet->step_ + verlet->step_rst_) % verlet->mdp.rstfreq == 0)
        {
            verlet->ucell.update_vel(verlet->vel);
            std::stringstream file;
            file << GlobalV::global_out_dir << "STRU_MD_" << verlet->step_ + verlet->step_rst_;
#ifdef __LCAO
            verlet->ucell.print_stru_file(GlobalC::ORB, file.str(), 1, 1);
#else
            verlet->ucell.print_stru_file(file.str(), 1, 1);
#endif
            verlet->write_restart();
        }

        verlet->step_++;
    }

    if (GlobalC::pot.out_potential == 2)
    {
        std::stringstream ssp;
        std::stringstream ssp_ave;
        ssp << GlobalV::global_out_dir << "ElecStaticPot";
        ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
        GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
    }

    ModuleBase::timer::tick("Run_MD_PW", "md_ions_pw");
    return;
}

void Run_MD_PW::md_cells_pw()
{
    ModuleBase::TITLE("Run_MD_PW", "md_cells_pw");
    ModuleBase::timer::tick("Run_MD_PW", "md_cells_pw");

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
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    GlobalC::pot.init_pot(0, GlobalC::pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    GlobalC::pot.newd();

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    //==================================================
    // create GlobalC::ppcell.tab_at , for trial wave functions.
    //==================================================
    GlobalC::wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if (GlobalV::NBANDS != 0 ) // liuyu update 2021-12-10
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
        throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
    }
#endif

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    // ion optimization begins
    // electron density optimization is included in ion optimization

    this->md_ions_pw();

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    ModuleBase::timer::tick("Run_MD_PW", "md_cells_pw");
}

void Run_MD_PW::md_force_virial(
    const int &istep,
    const int& numIon, 
    double &potential, 
    ModuleBase::Vector3<double>* force, 
    ModuleBase::matrix& virial)
{
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
        Electrons elec;
#ifdef __LCAO
        if (Exx_Global::Hybrid_Type::No == GlobalC::exx_global.info.hybrid_type)
        {
#endif
            elec.self_consistent(istep);
            eiter = elec.iter;
#ifdef __LCAO
        }
        else if (Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type)
        {
            throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        }
        else // Peize Lin add 2019-03-09
        {
            if (GlobalC::exx_global.info.separate_loop)
            {
                for (size_t hybrid_step = 0; hybrid_step != GlobalC::exx_global.info.hybrid_step; ++hybrid_step)
                {
                    elec.self_consistent(istep);
                    eiter += elec.iter;
                    if (elec.iter == 1 || hybrid_step == GlobalC::exx_global.info.hybrid_step - 1) // exx converge
                        break;
                    GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
                    GlobalC::exx_lip.cal_exx();
                }
            }
            else
            {
                elec.self_consistent(istep);
                eiter += elec.iter;
                GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
                elec.self_consistent(istep);
                eiter += elec.iter;
            }
        }
#endif
    }
    // mohan added 2021-01-28, perform stochastic calculations
    else if (GlobalV::CALCULATION == "md-sto")
    {
        Stochastic_Elec elec_sto;
        elec_sto.scf_stochastic(istep);
        eiter = elec_sto.iter;
    }

    ModuleBase::matrix fcs;
	Forces ff;
	ff.init(fcs);

	for(int ion=0;ion<numIon;ion++)
    {
		force[ion].x =fcs(ion, 0)/2.0;
		force[ion].y =fcs(ion, 1)/2.0;
		force[ion].z =fcs(ion, 2)/2.0;
	}

	if(GlobalV::STRESS)
	{
		Stress_PW ss;
		ss.cal_stress(virial);
        virial = 0.5 * virial;
	}

    potential = GlobalC::en.etot/2;
}