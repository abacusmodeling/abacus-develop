#include "run_md.h"

#include <vector>

#include "module_parameter/parameter.h"
#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "module_base/global_file.h"
#include "module_base/timer.h"
#include "module_io/print_info.h"
#include "module_plumed/plumed_interface.h"
#include "msst.h"
#include "nhchain.h"
#include "verlet.h"

namespace Run_MD
{

namespace
{
/// @brief collect the absolute Cartesian positions (in Bohr) of all the atoms
void collect_atomic_positions(const UnitCell& ucell, ModuleBase::Vector3<double>* positions)
{
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            /// the UnitCell stores tau in units of the lattice constant lat0
            positions[iat] = ucell.atoms[it].tau[ia] * ucell.lat0;
            ++iat;
        }
    }
}

/// @brief collect the cell matrix (in Bohr) in the row-major layout expected
/// by PLUMED (row i is the i-th lattice vector)
void collect_cell(const UnitCell& ucell, double* cell)
{
    /// the UnitCell stores the lattice vectors without the lattice constant
    const double* latvec = &ucell.latvec.e11;
    for (int i = 0; i < 9; ++i)
    {
        cell[i] = latvec[i] * ucell.lat0;
    }
}
} // namespace

void md_line(UnitCell& unit_in, ModuleESolver::ESolver* p_esolver, const Parameter& param_in)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::tick("Run_MD", "md_line");

    /// determine the md_type
    MD_base* mdrun;
    if (param_in.mdp.md_type == "fire")
    {
        mdrun = new FIRE(param_in, unit_in);
    }
    else if ((param_in.mdp.md_type == "nvt" && param_in.mdp.md_thermostat == "nhc") || param_in.mdp.md_type == "npt")
    {
        mdrun = new Nose_Hoover(param_in, unit_in);
    }
    else if (param_in.mdp.md_type == "nve" || param_in.mdp.md_type == "nvt")
    {
        mdrun = new Verlet(param_in, unit_in);
    }
    else if (param_in.mdp.md_type == "langevin")
    {
        mdrun = new Langevin(param_in, unit_in);
    }
    else if (param_in.mdp.md_type == "msst")
    {
        mdrun = new MSST(param_in, unit_in);
    }
    else
    {
        ModuleBase::WARNING_QUIT("md_line", "no such md_type!");
    }

    /// interface to PLUMED: no-op unless `plumed` is set in INPUT and ABACUS
    /// was compiled with PLUMED support (-DPLUMED_ROOT=...)
    ModulePlumed::PlumedInterface plumed;

    /// absolute Cartesian positions (Bohr) and cell (Bohr) handed over to
    /// PLUMED; note that the UnitCell stores tau and latvec in units of the
    /// lattice constant lat0 and that mdrun->pos only stores the displacement
    /// of the current MD step
    std::vector<ModuleBase::Vector3<double>> plumed_pos(unit_in.nat);
    double plumed_cell[9];

    /// md cycle
    while ((mdrun->step_ + mdrun->step_rst_) <= param_in.mdp.md_nstep && !mdrun->stop)
    {
        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, PARAM.globalv.global_readin_dir);

            if (param_in.mdp.plumed)
            {
                plumed.init(param_in.mdp.plumed_file,
                            unit_in.nat,
                            mdrun->get_md_dt(),
                            mdrun->allmass,
                            param_in.globalv.myrank);
                /// the initial configuration already carries positions and
                /// forces, so PLUMED can be called right away
                collect_atomic_positions(unit_in, plumed_pos.data());
                collect_cell(unit_in, plumed_cell);
                plumed.move(mdrun->step_ + mdrun->step_rst_,
                            plumed_pos.data(),
                            mdrun->force,
                            plumed_cell,
                            mdrun->potential,
                            param_in.inp.cal_stress ? mdrun->virial.c : nullptr);
            }
        }
        else
        {
            ModuleIO::print_screen(0, 0, mdrun->step_ + mdrun->step_rst_);
            mdrun->first_half(GlobalV::ofs_running);

            /// update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver,
                                  mdrun->step_,
                                  unit_in,
                                  mdrun->potential,
                                  mdrun->force,
                                  param_in.inp.cal_stress,
                                  mdrun->virial);

            if (plumed.active())
            {
                /// let PLUMED compute the collective variables and add the
                /// biasing forces to mdrun->force; they are used by the second
                /// half of the current step and by the first half of the next
                collect_atomic_positions(unit_in, plumed_pos.data());
                collect_cell(unit_in, plumed_cell);
                plumed.move(mdrun->step_ + mdrun->step_rst_,
                            plumed_pos.data(),
                            mdrun->force,
                            plumed_cell,
                            mdrun->potential,
                            param_in.inp.cal_stress ? mdrun->virial.c : nullptr);
            }

            mdrun->second_half();

            MD_func::compute_stress(unit_in,
                                    mdrun->vel,
                                    mdrun->allmass,
                                    param_in.inp.cal_stress,
                                    mdrun->virial,
                                    mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     unit_in.nat,
                                                     mdrun->frozen_freedom_,
                                                     mdrun->allmass,
                                                     mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_dumpfreq == 0)
        {
            mdrun->print_md(GlobalV::ofs_running, PARAM.inp.cal_stress);

            MD_func::dump_info(mdrun->step_ + mdrun->step_rst_,
                               PARAM.globalv.global_out_dir,
                               unit_in,
                               param_in,
                               mdrun->virial,
                               mdrun->force,
                               mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_restartfreq == 0)
        {
            unit_in.update_vel(mdrun->vel);
            std::stringstream file;
            file << PARAM.globalv.global_stru_dir << "STRU_MD_" << mdrun->step_ + mdrun->step_rst_;
            // changelog 20240509
            // because I move out the dependence on GlobalV from UnitCell::print_stru_file
            // so its parameter is calculated here
            bool need_orb = PARAM.inp.basis_type=="pw";
            need_orb = need_orb && PARAM.inp.init_wfc.substr(0, 3)=="nao";
            need_orb = need_orb || PARAM.inp.basis_type=="lcao";
            need_orb = need_orb || PARAM.inp.basis_type=="lcao_in_pw";
            unit_in.print_stru_file(file.str(), 
                                    PARAM.inp.nspin, 
                                    false, // Cartesian coordinates
                                    PARAM.inp.calculation == "md", 
                                    PARAM.inp.out_mul,
                                    need_orb,
                                    PARAM.globalv.deepks_setorb,
                                    GlobalV::MY_RANK);
            mdrun->write_restart(PARAM.globalv.global_out_dir);
        }

        mdrun->step_++;
    }

    /// flush the COLVAR and other PLUMED output files
    plumed.finalize();

    delete mdrun;

    ModuleBase::timer::tick("Run_MD", "md_line");
    return;
}

} // namespace Run_MD