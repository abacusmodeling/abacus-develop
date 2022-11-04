#include "mdrun.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"
#include "module_esolver/esolver.h"
#include "../src_io/print_info.h"

MDrun::MDrun(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in):
    mdp(MD_para_in),
    ucell(unit_in)
{
    if(mdp.md_seed >= 0)
    {
        srand(mdp.md_seed);
    }

    stop = false;

    allmass = new double [ucell.nat];
    pos = new ModuleBase::Vector3<double> [ucell.nat];
    vel = new ModuleBase::Vector3<double> [ucell.nat];
    ionmbl = new ModuleBase::Vector3<int> [ucell.nat];
    force = new ModuleBase::Vector3<double> [ucell.nat];
    virial.create(3,3);
    stress.create(3,3);

    // convert to a.u. unit
    mdp.md_dt /= ModuleBase::AU_to_FS;
    mdp.md_tfirst /= ModuleBase::Hartree_to_K;
    mdp.md_tlast /= ModuleBase::Hartree_to_K;
    mdp.md_tfreq *= ModuleBase::AU_to_FS;

    step_ = 0;
    step_rst_ = 0;

    MD_func::InitPos(ucell, pos);
    MD_func::InitVel(ucell, mdp.md_tfirst, allmass, frozen_freedom_, ionmbl, vel);
}

MDrun::~MDrun()
{
    delete []allmass;
    delete []pos;
    delete []vel;
    delete []ionmbl;
    delete []force;
}

void MDrun::setup(ModuleESolver::ESolver *p_esolver)
{
    if(mdp.md_restart)
    {
        restart();
    }

    Print_Info::print_screen(0, 0, step_ + step_rst_);

    MD_func::force_virial(p_esolver, step_, mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    temperature_ = 2*kinetic/(double(3*ucell.nat-frozen_freedom_))*ModuleBase::Hartree_to_K;
}

void MDrun::first_half()
{
    if(GlobalV::MY_RANK==0) //only first rank do md
    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                vel[i][k] += 0.5*force[i][k]*mdp.md_dt/allmass[i];
                pos[i][k] += vel[i][k]*mdp.md_dt;
            }
        }
    }
#ifdef __MPI
    MPI_Bcast(pos , ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(vel , ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();
    MD_func::InitPos(ucell, pos);
}

void MDrun::second_half()
{
    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                vel[i][k] += 0.5*force[i][k]*mdp.md_dt/allmass[i];
            }
        }
    }
}

void MDrun::outputMD(std::ofstream &ofs, bool cal_stress)
{
    if(GlobalV::MY_RANK) return;

    temperature_ = 2*kinetic/(double(3*ucell.nat-frozen_freedom_))*ModuleBase::Hartree_to_K;

    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
    double press = 0.0;
    for(int i=0;i<3;i++)
    {
        press += stress(i,i)/3;
    }

    std::cout << " ------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << " " << std::left << std::setw(20) << "Energy" 
            << std::left << std::setw(20) << "Potential" 
            << std::left << std::setw(20) << "Kinetic" 
            << std::left << std::setw(20) << "Temperature";
    if(cal_stress)
    {
        std::cout << std::left << std::setw(20) << "Pressure (KBAR)";
    }
    std::cout << std::endl;
    std::cout << " " << std::left << std::setw(20) << potential+kinetic
            << std::left << std::setw(20) << potential
            << std::left << std::setw(20) << kinetic
            << std::left << std::setw(20) << temperature_;
    if(cal_stress)
    {
        std::cout << std::left << std::setw(20) << press*unit_transform;
    }
    std::cout << std::endl;
	std::cout << " ------------------------------------------------------------------------------------------------" << std::endl;

    ofs.unsetf(ios::fixed);
    ofs << std::setprecision(8) << std::endl;
    ofs << std::endl;
    ofs << " ------------------------------------------------------------------------------------------------" << std::endl;
	ofs << " " << std::left << std::setw(20) << "Energy" 
        << std::left << std::setw(20) << "Potential" 
        << std::left << std::setw(20) << "Kinetic" 
        << std::left << std::setw(20) << "Temperature"; 
    if(cal_stress)
    {
        ofs << std::left << std::setw(20) << "Pressure (KBAR)";
    }
    ofs << std::endl;
    ofs << " " << std::left << std::setw(20) << potential+kinetic
        << std::left << std::setw(20) << potential
        << std::left << std::setw(20) << kinetic
        << std::left << std::setw(20) << temperature_;
    if(cal_stress)
    {
        ofs << std::left << std::setw(20) << press*unit_transform;
    }
    ofs << std::endl;
    ofs << " ------------------------------------------------------------------------------------------------" << std::endl;
    if(cal_stress)
    {
        MD_func::outStress(virial, stress);
    }
    ofs << std::endl;
    ofs << std::endl;
}

void MDrun::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
		file.close();
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void MDrun::restart()
{
    bool ok = true;

    if(!GlobalV::MY_RANK)
    {
        std::stringstream ssc;
        ssc << GlobalV::global_readin_dir << "Restart_md.dat";
        std::ifstream file(ssc.str().c_str());

        if(!file)
        {
            ok = false;
        }

        if(ok)
        {
            file >> step_rst_;
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if(!ok)
    {
        ModuleBase::WARNING_QUIT("mdrun", "no Restart_md.dat !");
    }

#ifdef __MPI
    MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}