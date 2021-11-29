#include "NVE.h"
#include "MD_func.h"

NVE::NVE(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{}

NVE::~NVE(){}

void NVE::setup()
{
    ModuleBase::TITLE("NVE", "setup");
    ModuleBase::timer::tick("NVE", "setup");

    MD_func::force_virial(mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    ModuleBase::timer::tick("NVE", "setup");
}

void NVE::first_half()
{
    ModuleBase::TITLE("NVE", "first_half");
    ModuleBase::timer::tick("NVE", "first_half");

    for(int i=0; i<ucell.nat; ++i)
    {
        if(ionmbl[i].x)
        {
            pos[i].x += vel[i].x*mdp.dt + 0.5*force[i].x*mdp.dt*mdp.dt/allmass[i];
            vel[i].x += 0.5*force[i].x*mdp.dt/allmass[i];
        }
        if(ionmbl[i].y)
        {
            pos[i].y += vel[i].y*mdp.dt + 0.5*force[i].y*mdp.dt*mdp.dt/allmass[i];
            vel[i].y += 0.5*force[i].y*mdp.dt/allmass[i];
        }
        if(ionmbl[i].z)
        {
            pos[i].z += vel[i].z*mdp.dt + 0.5*force[i].z*mdp.dt*mdp.dt/allmass[i];
            vel[i].z += 0.5*force[i].z*mdp.dt/allmass[i];
        }
    }

    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();

    ModuleBase::timer::tick("NVE", "first_half");
}

void NVE::second_half()
{
    ModuleBase::TITLE("NVE", "second_half");
    ModuleBase::timer::tick("NVE", "second_half");

    for(int i=0; i<ucell.nat; ++i)
    {
        if(ionmbl[i].x)
        {
            vel[i].x += 0.5*force[i].x*mdp.dt/allmass[i];
        }
        if(ionmbl[i].y)
        {
            vel[i].y += 0.5*force[i].y*mdp.dt/allmass[i];
        }
        if(ionmbl[i].z)
        {
            vel[i].z += 0.5*force[i].z*mdp.dt/allmass[i];
        }
    }

    ModuleBase::timer::tick("NVE", "second_half");
}

void NVE::outputMD()
{
    ModuleBase::TITLE("NVE", "outputMD");
    ModuleBase::timer::tick("NVE", "outputMD");

    temperature_ = 2*kinetic/(double(3*ucell.nat-frozen_freedom_))*ModuleBase::Hartree_to_K;

    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
    double press = 0.0;
    for(int i=0;i<3;i++)
    {
        press += stress(i,i)/3;
    }

    std::cout << std::setprecision(6) << std::setiosflags(ios::showpos) << std::setiosflags(ios::fixed);
    std::cout << " " << std::left << std::setw(20) << "Energy" 
            << std::left << std::setw(20) << "Potential" 
            << std::left << std::setw(20) << "Kinetic" 
            << std::left << std::setw(20) << "Temperature" 
            << std::left << std::setw(20) << "Pressure (KBAR)" <<std::endl;
    std::cout << " " << std::left << std::setw(20) << potential+kinetic
            << std::left << std::setw(20) << potential
            << std::left << std::setw(20) << kinetic
            << std::left << std::setw(20) << temperature_
            << std::left << std::setw(20) << press*unit_transform <<std::endl;
	std::cout << " ------------------------------------------------------------" << std::endl;

    MD_func::outStress(virial, stress);

    ModuleBase::timer::tick("NVE", "outputMD");
}

void NVE::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ << std::endl;
		file.close();
	}
}

void NVE::restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ifstream file(ssc.str().c_str());

        if(!file)
		{
			std::cout<< "please ensure whether 'Restart_md.dat' exists!" << std::endl;
            ModuleBase::WARNING_QUIT("MSST", "no Restart_md.dat ！");
		}

		file >> step_rst_;

		file.close();
	}

#ifdef __MPI
	MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    step_ = step_rst_;
}