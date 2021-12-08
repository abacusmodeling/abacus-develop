#include "NVT_NHC.h"
#include "MD_func.h"

NVT_NHC::NVT_NHC(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    // convert to a.u. unit
    mdp.tfreq *= ModuleBase::AU_to_FS;

    // init NHC
    Q = new double [mdp.MNHC];
	G = new double [mdp.MNHC];
	eta = new double [mdp.MNHC];
	veta = new double [mdp.MNHC];

    for(int i=0; i<mdp.MNHC; ++i)
    {
        eta[i] = veta[i] = G[i] = 0;
    }

    //w[0] = 1;

    w[0] = 0.784513610477560;
	w[6] = 0.784513610477560;
	w[1] = 0.235573213359357;
	w[5] = 0.235573213359357;
	w[2] = -1.17767998417887;
	w[4] = -1.17767998417887;
	w[3] = 1-w[0]-w[1]-w[2]-w[4]-w[5]-w[6];
}

NVT_NHC::~NVT_NHC()
{
    delete []Q;
    delete []G;
    delete []eta;
    delete []veta;
}

void NVT_NHC::setup()
{
    ModuleBase::TITLE("NVT_NHC", "setup");
    ModuleBase::timer::tick("NVT_NHC", "setup");

    Verlet::setup();

    temp_target();
    
    update_mass();
    
    for(int m=1; m<mdp.MNHC; ++m)
    {
        G[m] = (Q[m-1]*veta[m-1]*veta[m-1]-t_target) / Q[m];
    }

    ModuleBase::timer::tick("NVT_NHC", "setup");
}

void NVT_NHC::first_half()
{
    ModuleBase::TITLE("NVT_NHC", "first_half");
    ModuleBase::timer::tick("NVT_NHC", "first_half");

    // update target T
    temp_target();

    integrate();

    Verlet::first_half();

    ModuleBase::timer::tick("NVT_NHC", "first_half");
}

void NVT_NHC::second_half()
{
    ModuleBase::TITLE("NVT_NHC", "second_half");
    ModuleBase::timer::tick("NVT_NHC", "second_half");

    Verlet::second_half();

    integrate();

    ModuleBase::timer::tick("NVT_NHC", "second_half");
}

void NVT_NHC::outputMD()
{
    Verlet::outputMD();
}

void NVT_NHC::write_restart()
{
    Verlet::write_restart();
}

void NVT_NHC::restart()
{
    std::cout << "NHC" << std::endl;
    Verlet::restart();
}

void NVT_NHC::integrate()
{
    double scale = 1.0;
    kinetic = MD_func::GetAtomKE(ucell.nat, vel, allmass);
    double KE = kinetic;

    // update mass
    update_mass();

    // update force
    if(Q[0] > 0) 
    {
        G[0] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_target) / Q[0];
    }
    else 
    {
        G[0] = 0;
    }

    for(int i=0; i<nc; ++i)
    {
        for(int j=0; j<nsy; ++j)
        {
            double delta = w[j] * mdp.dt / nc;

            // propogate veta
            veta[mdp.MNHC-1] += G[mdp.MNHC-1] * delta /4.0;

            for(int m=mdp.MNHC-2; m>=0; --m)
            {
                double aa = exp(-veta[m]*delta/8.0);
                veta[m] = veta[m] * aa * aa + G[m] * aa * delta /4.0;
            }

            scale *= exp(-veta[0]*delta/2.0);
            KE = kinetic * scale * scale;

            // update force
            if(Q[0] > 0) 
            {
                G[0] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_target) / Q[0];
            }
            else 
            {
                G[0] = 0;
            }

            // propogate eta
            for(int m=0; m<mdp.MNHC; ++m)
            {
                eta[m] += veta[m] * delta / 2.0;
            }

            // propogate veta
            for(int m=0; m<mdp.MNHC-1; ++m)
            {
                double aa = exp(-veta[m+1]*delta/8.0);
                veta[m] = veta[m] * aa * aa + G[m] * aa * delta /4.0;

                G[m+1] = (Q[m]*veta[m]*veta[m]-t_target) / Q[m+1];
            }
            veta[mdp.MNHC-1] += G[mdp.MNHC-1] * delta /4.0;
        }
    }
    
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
}

void NVT_NHC::temp_target()
{
    double delta = (double)step_ / GlobalV::NSTEP;
    t_target = mdp.tfirst + delta * (mdp.tlast - mdp.tfirst);
}

void NVT_NHC::update_mass()
{
    Q[0] = (3*ucell.nat - frozen_freedom_) * t_target / mdp.tfreq / mdp.tfreq;
    for(int m=1; m<mdp.MNHC; ++m)
    {
        Q[m] = t_target / mdp.tfreq / mdp.tfreq;
    }
}