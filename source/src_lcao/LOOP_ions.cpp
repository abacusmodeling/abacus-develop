#include "LOOP_ions.h"
#include "../src_pw/global.h"
#include "../module_orbital/parallel_orbitals.h"
#include "../src_pdiag/pdiag_double.h"
#include "../module_base/global_function.h"
#include "../src_io/write_HS.h"
#include "../src_io/print_info.h"
#include "../src_io/cal_r_overlap_R.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"
#include "dmft.h"
#include "src_lcao/LCAO_matrix.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"    //caoyu add 2021-07-26
#endif

LOOP_ions::LOOP_ions()
{}

LOOP_ions::~LOOP_ions()
{}

void LOOP_ions::opt_ions(ModuleESolver::ESolver* p_esolver)
{
    ModuleBase::TITLE("LOOP_ions", "opt_ions");
    ModuleBase::timer::tick("LOOP_ions", "opt_ions");

	if(GlobalV::OUT_LEVEL=="i")
	{
		std::cout << std::setprecision(12);
    	std::cout<< " " << std::setw(7)<< "ISTEP" 
		<<std::setw(5)<< "NE"
		<<std::setw(15)<< "ETOT(eV)"
		<<std::setw(15)<< "EDIFF(eV)"
        <<std::setw(15)<< "MAX_F(eV/A)"
        <<std::setw(15)<< "TRADIUS(Bohr)"
		<<std::setw(8)<< "UPDATE"
		<<std::setw(11)<< "ETIME(MIN)"
		<<std::setw(11)<< "FTIME(MIN)"
        <<std::endl;
	}

    // Geometry optimization algorithm setup.
    if (GlobalV::CALCULATION=="relax")
    {
        //Ions_Move_Methods
        IMM.allocate();
    }
    if (GlobalV::CALCULATION=="cell-relax")
    {
        //Ions_Move_Methods
        IMM.allocate();
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    } 

    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;
    while (istep <= GlobalV::RELAX_NMAX && !stop)
    {
        time_t estart = time(NULL);

        // xiaohui add "m" option, 2015-09-16
        if (GlobalV::OUT_LEVEL == "ie" || GlobalV::OUT_LEVEL == "m")
        {
            Print_Info::print_screen(stress_step, force_step, istep);
        }

        // solve electronic structures in terms of LCAO
        p_esolver->Run(this->istep - 1, GlobalC::ucell);

        time_t eend = time(NULL);

        //for second-order extrapolation
        if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            CE.update_all_pos(GlobalC::ucell);
        }

        time_t fstart = time(NULL);
        if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            stop = this->force_stress(istep, force_step, stress_step, p_esolver);
        }
        time_t fend = time(NULL);

        // PLEASE move the details of CE to other places
        // mohan add 2021-03-25
        //xiaohui add 2014-07-07, for second-order extrapolation
        if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            CE.save_pos_next(GlobalC::ucell);
        }

		if(GlobalV::OUT_LEVEL=="i")
		{
			double etime_min = difftime(eend, estart)/60.0; 
			double ftime_min = difftime(fend, fstart)/60.0; 
			std::stringstream ss;
			ss << GlobalV::RELAX_METHOD << istep;
			
			std::cout << " " << std::setw(7) << ss.str() 
			<< std::setw(5) << p_esolver->getniter()
			<< std::setw(15) << std::setprecision(6) << GlobalC::en.etot * ModuleBase::Ry_to_eV 
			<< std::setw(15) << IMM.get_ediff() * ModuleBase::Ry_to_eV
			<< std::setprecision(3)
			<< std::setw(15) << IMM.get_largest_grad() * ModuleBase::Ry_to_eV / 0.529177
			<< std::setw(15) << IMM.get_trust_radius()
			<< std::setw(8) << IMM.get_update_iter()
			<< std::setprecision(2) << std::setw(11) << etime_min
			<< std::setw(11) << ftime_min << std::endl;
		}

        ++istep;
    }

    ModuleBase::timer::tick("LOOP_ions", "opt_ions");
    return;
}


bool LOOP_ions::force_stress(
    const int& istep,
    int& force_step,
    int& stress_step,
    ModuleESolver::ESolver* p_esolver)
{
    ModuleBase::TITLE("LOOP_ions", "force_stress");

    if (!GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {
        return 1;
    }
    ModuleBase::timer::tick("LOOP_ions", "force_stress");

    // set force matrix
    ModuleBase::matrix fcs;
    // set stress matrix
    ModuleBase::matrix scs;

    p_esolver->cal_Force(fcs);
    p_esolver->cal_Stress(scs);

    //--------------------------------------------------
    // only forces are needed, no stresses are needed
    //--------------------------------------------------
    if (GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {

#ifdef __MPI
        atom_arrange::delete_vector(
            GlobalV::ofs_running,
            GlobalV::SEARCH_PBC,
            GlobalC::GridD,
            GlobalC::ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);
#endif

        if (GlobalV::CALCULATION == "relax")
        {
            IMM.cal_movement(istep, istep, fcs, GlobalC::en.etot);

            if (IMM.get_converged() || (istep == GlobalV::RELAX_NMAX))
            {
                ModuleBase::timer::tick("LOOP_ions", "force_stress");
                return 1; // 1 means converged
            }
            else // ions are not converged
            {
                CE.update_istep(istep);
                CE.extrapolate_charge();

                if (GlobalC::pot.chg_extrap == "dm")
                {
                }
                else
                {
                    GlobalC::pot.init_pot(istep, GlobalC::sf.strucFac);
                }
            }
            ModuleBase::timer::tick("LOOP_ions", "force_stress");
            return 0;
        }
        else
        {
            ModuleBase::timer::tick("LOOP_ions", "force_stress");
            return 1;
        }

    }

    //    static bool converged_force = false;
    static bool converged_stress = false;

    if (!GlobalV::CAL_FORCE && GlobalV::CAL_STRESS)
    {

#ifdef __MPI
        atom_arrange::delete_vector(
            GlobalV::ofs_running,
            GlobalV::SEARCH_PBC,
            GlobalC::GridD,
            GlobalC::ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);
#endif
        if (GlobalV::CALCULATION == "cell-relax")
        {
            LCM.cal_lattice_change(stress_step, scs, GlobalC::en.etot);
            converged_stress = LCM.get_converged();
            if (converged_stress)
            {
                ModuleBase::timer::tick("LOOP_ions", "force_stress");
                return 1;
            }
            else
            {
                Variable_Cell::init_after_vc(p_esolver);
                GlobalC::pot.init_pot(stress_step, GlobalC::sf.strucFac);

                ++stress_step;
                ModuleBase::timer::tick("LOOP_ions", "force_stress");
                return 0;
            }
        }
        else
        {
            ModuleBase::timer::tick("LOOP_ions", "force_stress");
            return 1;
        }
    }

    if (GlobalV::CAL_FORCE && GlobalV::CAL_STRESS)
    {
        atom_arrange::delete_vector(
            GlobalV::ofs_running,
            GlobalV::SEARCH_PBC,
            GlobalC::GridD,
            GlobalC::ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);

        if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            IMM.cal_movement(istep, force_step, fcs, GlobalC::en.etot);

            if (IMM.get_converged())
            {
                force_step = 1;


                if (GlobalV::CALCULATION == "cell-relax")
                {
                    LCM.cal_lattice_change(stress_step, scs, GlobalC::en.etot);
                    converged_stress = LCM.get_converged();
                    if (converged_stress)
                    {
                        ModuleBase::timer::tick("LOOP_ions", "force_stress");
                        return 1;
                    }
                    else
                    {
                        Variable_Cell::init_after_vc(p_esolver);
                        GlobalC::pot.init_pot(stress_step, GlobalC::sf.strucFac);

                        ++stress_step;
                        ModuleBase::timer::tick("LOOP_ions", "force_stress");
                        return 0;
                    }
                }
                else
                {
                    ModuleBase::timer::tick("LOOP_ions", "force_stress");
                    return 1;
                }

            }
            else
            {
                CE.update_istep(force_step);
                CE.extrapolate_charge();

                if (GlobalC::pot.chg_extrap == "dm")
                {
                }
                else
                {
                    GlobalC::pot.init_pot(istep, GlobalC::sf.strucFac);
                }
                ++force_step;
                ModuleBase::timer::tick("LOOP_ions", "force_stress");
                return 0;
            }
        }
        else
        {
            ModuleBase::timer::tick("LOOP_ions", "force_stress");
            return 1;
        }
    }
    ModuleBase::timer::tick("LOOP_ions", "force_stress");
    return 0;
}


