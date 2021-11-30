#include "NVT_ADS.h"
#include "MD_func.h"

NVT_ADS::NVT_ADS(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    // convert to a.u. unit
    mdp.viscosity *= ModuleBase::AU_to_FS;

    nraise = mdp.viscosity * mdp.dt;
}

NVT_ADS::~NVT_ADS(){}

void NVT_ADS::setup()
{
    ModuleBase::TITLE("NVT_ADS", "setup");
    ModuleBase::timer::tick("NVT_ADS", "setup");

    if(mdp.rstMD)
    {
        restart();
    }

    MD_func::force_virial(mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    ModuleBase::timer::tick("NVT_ADS", "setup");
}

void NVT_ADS::first_half()
{
    ModuleBase::TITLE("NVT_ADS", "first_half");
    ModuleBase::timer::tick("NVT_ADS", "first_half");

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

    ModuleBase::timer::tick("NVT_ADS", "first_half");
}

void NVT_ADS::second_half()
{
    ModuleBase::TITLE("NVT_ADS", "second_half");
    ModuleBase::timer::tick("NVT_ADS", "second_half");

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

    ModuleBase::timer::tick("NVT_ADS", "second_half");
}

void NVT_ADS::outputMD()
{
    Verlet::outputMD();
}

void NVT_ADS::write_restart()
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

void NVT_ADS::restart()
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