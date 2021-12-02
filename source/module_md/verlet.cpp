#include "verlet.h"
#include "MD_func.h"

Verlet::Verlet(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in):
    mdp(MD_para_in),
    ucell(unit_in)
{
    allmass = new double [ucell.nat];
    pos = new ModuleBase::Vector3<double> [ucell.nat];
    vel = new ModuleBase::Vector3<double> [ucell.nat];
    ionmbl = new ModuleBase::Vector3<int> [ucell.nat];
    force = new ModuleBase::Vector3<double> [ucell.nat];
    virial.create(3,3);
    stress.create(3,3);

    // convert to a.u. unit
    mdp.dt/=ModuleBase::AU_to_FS;
	temperature_=mdp.tfirst/ModuleBase::Hartree_to_K;

    // LJ parameters
    mdp.rcut_lj *= ModuleBase::ANGSTROM_AU;
	mdp.epsilon_lj /= ModuleBase::Hartree_to_eV;
	mdp.sigma_lj *= ModuleBase::ANGSTROM_AU;

    step_ = 0;
    step_rst_ = 0;
    if(mdp.rstMD) unit_in.set_vel = 1;

    //frozen_freedom_ = MD_func::getMassMbl(ucell, mdp, allmass, ionmbl);
    MD_func::InitPos(ucell, pos);
    MD_func::InitVel(ucell, temperature_, allmass, frozen_freedom_, ionmbl, vel);
}

Verlet::~Verlet()
{
    delete []allmass;
    delete []pos;
    delete []vel;
    delete []ionmbl;
    delete []force;
}

void Verlet::setup()
{
    if(mdp.rstMD)
    {
        restart();
    }

    MD_func::force_virial(mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;
}

void Verlet::first_half()
{
    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                pos[i][k] += vel[i][k]*mdp.dt + 0.5*force[i][k]*mdp.dt*mdp.dt/allmass[i];
                vel[i][k] += 0.5*force[i][k]*mdp.dt/allmass[i];
            }
        }
    }

    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();
}

void Verlet::second_half()
{
    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                vel[i][k] += 0.5*force[i][k]*mdp.dt/allmass[i];
            }
        }
    }
}

void Verlet::outputMD()
{
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
}

void Verlet::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ << std::endl;
		file.close();
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Verlet::restart()
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