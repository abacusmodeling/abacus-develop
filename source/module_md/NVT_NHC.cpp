#include "NVT_NHC.h"
#include "MD_func.h"

NVT_NHC::NVT_NHC(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    // convert to a.u. unit
    if(mdp.Qmass < ucell.nat) mdp.Qmass = ucell.nat;
    mdp.Qmass /= ModuleBase::AU_to_MASS;

    // init NHC
    Q = new double [mdp.MNHC];
	G = new double [mdp.MNHC];
	NHCeta = new double [mdp.MNHC];
	NHCveta = new double [mdp.MNHC];

    Q[0] = mdp.Qmass;
    for(int i=0; i<mdp.MNHC; ++i)
    {
        if(i>0) Q[i] = mdp.Qmass / (3*ucell.nat - frozen_freedom_);
        NHCeta[i] = rand()/double(RAND_MAX) - 0.5;
        NHCveta[i] = MD_func::gaussrand() * sqrt(t_last / Q[i]);
    }

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
    delete []NHCeta;
    delete []NHCveta;
}

void NVT_NHC::setup()
{
    ModuleBase::TITLE("NVT_NHC", "setup");
    ModuleBase::timer::tick("NVT_NHC", "setup");

    Verlet::setup();

    G[0] = (2*kinetic - (3*ucell.nat - frozen_freedom_)*t_last) / Q[0];
    for(int m=1; m<mdp.MNHC; ++m)
    {
        G[m] = (Q[m-1]*NHCveta[m-1]*NHCveta[m-1]-t_last) / Q[m];
    }

    ModuleBase::timer::tick("NVT_NHC", "setup");
}

void NVT_NHC::first_half()
{
    ModuleBase::TITLE("NVT_NHC", "first_half");
    ModuleBase::timer::tick("NVT_NHC", "first_half");

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

    // update force
    G[0] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_last) / Q[0];
    // for(int m=1; m<mdp.MNHC; ++m)
    // {
    //     G[m] = (Q[m-1]*NHCveta[m-1]*NHCveta[m-1]-t_last) / Q[m];
    // }

    for(int i=0; i<nc; ++i)
    {
        for(int j=0; j<nsy; ++j)
        {
            double delta = w[j] * mdp.dt / nc;

            // propogate veta
            NHCveta[mdp.MNHC-1] += G[mdp.MNHC-1] * delta /4.0;

            for(int m=mdp.MNHC-2; m>=0; --m)
            {
                double aa = exp(-NHCveta[m]*delta/8.0);
                NHCveta[m] = NHCveta[m] * aa * aa + G[m] * aa * delta /4.0;
            }

            scale *= exp(-NHCveta[0]*delta/2.0);
            KE = kinetic * scale * scale;

            // propogate eta
            for(int m=0; m<mdp.MNHC; ++m)
            {
                NHCeta[m] += NHCveta[m] * delta / 2.0;
            }

            // propogate veta
            for(int m=0; m<mdp.MNHC-1; ++m)
            {
                if(m==0)
                {
                    G[m] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_last) / Q[m];
                }
                else
                {
                    G[m] = (Q[m-1]*NHCveta[m-1]*NHCveta[m-1]-t_last) / Q[m];
                }
                
                double aa = exp(-NHCveta[m+1]*delta/8.0);
                NHCveta[m] = NHCveta[m] * aa * aa + G[m] * aa * delta /4.0;
            }
            G[mdp.MNHC-1] = (Q[mdp.MNHC-2]*NHCveta[mdp.MNHC-2]*NHCveta[mdp.MNHC-2]-t_last) / Q[mdp.MNHC-1];
            NHCveta[mdp.MNHC-1] += G[mdp.MNHC-1] * delta /4.0;
        }
    }
    
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
}