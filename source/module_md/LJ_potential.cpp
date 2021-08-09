#include "LJ_potential.h"
#include "../input.h"

LJ_potential::LJ_potential(){}

LJ_potential::~LJ_potential(){}

double LJ_potential::Lennard_Jones(UnitCell_pseudo &ucell_c, 
                    Grid_Driver &grid_neigh, 
                    Vector3<double> *force, 
                    matrix &stress)
{
    double distance, potential = 0; //initialize
    int index = 0;
    Vector3<double> tau1, tau2, dtau;
    for(int it=0; it<ucell_c.ntype; ++it)
    {
        Atom* atom1 = &ucell_c.atoms[it];
        for(int ia=0; ia<atom1->na; ++ia)
	    {	
		    force[index].set(0,0,0);   //initialize
            tau1 = atom1->tau[ia];
            grid_neigh.Find_atom(ucell_c, tau1, it, ia);
            for(int ad=0; ad<grid_neigh.getAdjacentNum(); ++ad)
            {
                tau2 = grid_neigh.getAdjacentTau(ad);
                dtau = (tau1 - tau2) * ucell_c.lat0;
                distance = dtau.norm();
                if(distance <= INPUT.mdp.rcut_lj)
			    {
				    potential += LJ_energy(distance) - LJ_energy(INPUT.mdp.rcut_lj);
				    force[index] += LJ_force(distance, dtau);
			    }
            }
            index++;
	    }
    }
	return potential/2.0;
}

double LJ_potential::LJ_energy(const double d)
{
	return 4*INPUT.mdp.epsilon_lj*( pow(INPUT.mdp.sigma_lj/d, 12) - pow(INPUT.mdp.sigma_lj/d, 6) );
}

Vector3<double> LJ_potential::LJ_force(const double d, const Vector3<double> dr)
{
	double coff = 4*INPUT.mdp.epsilon_lj*( 12*pow(INPUT.mdp.sigma_lj/d, 12) - 6*pow(INPUT.mdp.sigma_lj/d, 6) )/pow(d,2);
	return dr*coff;
}