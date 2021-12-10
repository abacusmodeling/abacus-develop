#include "FIRE.h"
#include "MD_func.h"

FIRE::FIRE(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    dt_max = -1.0;
    alpha_start = 0.10;
    alpha = alpha_start;

    finc = 1.1;
    fdec = 0.5;
    f_alpha = 0.99;
    N_min = 4;
    negative_count = 0;
}

FIRE::~FIRE(){}

void FIRE::setup()
{
    ModuleBase::TITLE("FIRE", "setup");
    ModuleBase::timer::tick("FIRE", "setup");
    
    Verlet::setup();

    check_force();

    ModuleBase::timer::tick("FIRE", "setup");
}

void FIRE::first_half()
{
    ModuleBase::TITLE("FIRE", "first_half");
    ModuleBase::timer::tick("FIRE", "first_half");

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

    check_FIRE();

    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                pos[i][k] += vel[i][k]*mdp.dt;
            }
        }
    }

    check_force();

    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();
    MD_func::InitPos(ucell, pos);

    ModuleBase::timer::tick("FIRE", "first_half");
}

void FIRE::second_half()
{
    ModuleBase::TITLE("FIRE", "second_half");
    ModuleBase::timer::tick("FIRE", "second_half");

    Verlet::second_half();

    ModuleBase::timer::tick("FIRE", "second_half");
}

void FIRE::outputMD()
{
    Verlet::outputMD();

    GlobalV::ofs_running << " LARGEST GRAD (eV/A)  : " 
        << max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU << std::endl;
    std::cout << " LARGEST GRAD (eV/A)  : " 
        << max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU << std::endl;
}

void FIRE::write_restart()
{
    Verlet::write_restart();
}

void FIRE::restart()
{
    Verlet::restart();
}

void FIRE::check_force()
{
    max = 0;

    for(int i=0; i<ucell.nat; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            if( max < abs(force[i][j]) )
            {
                max = abs(force[i][j]);
            }
        }
    }

    if(2.0*max < GlobalV::FORCE_THR)
    {
        stop = true;
    }
}

void FIRE::check_FIRE()
{
	double P = 0;
	double sumforce = 0;
	double normvel = 0;

    if(dt_max < 0) dt_max = 2.5 * mdp.dt;   //initial dt_max
	
	for(int i =0; i<ucell.nat; ++i)
    {
		P += vel[i].x*force[i].x+vel[i].y*force[i].y+vel[i].z*force[i].z;
		sumforce += force[i].norm2();
		normvel += vel[i].norm2();
	}
	
	sumforce = sqrt(sumforce);
	normvel = sqrt(normvel);
	
	for(int i=0; i<ucell.nat; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            vel[i][j] = (1.0-alpha) * vel[i][j] + alpha * force[i][j] / sumforce * normvel;
        }
	}
	
	if(P > 0)
	{
	    negative_count++;
		if(negative_count >= N_min)
        {
			mdp.dt = min(mdp.dt*finc, dt_max);
			alpha *= f_alpha;
		}
	}
	else
	{
		mdp.dt *= fdec;
		negative_count = 0;
		
		for(int i=0; i<ucell.nat; ++i)
        {
			for(int j=0; j<3; ++j)
            {
                vel[i][j] = 0;
            }
		}
		
		alpha = alpha_start;
	}
}

