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

double LJ_potential::Lennard_Jones(UnitCell_pseudo &ucell_c, 
                    CMD_neighbor &cmd_neigh,
                    Vector3<double> *force, 
                    matrix &stress)
{
    double potential = 0; //initialize

	for(int i=0; i<ucell_c.nat; i++)
	{	
		force[i].set(0,0,0);   //initialize

        int T1 = ucell_c.iat2it[i];
        int I1 = ucell_c.iat2ia[i];
        Vector3<double> taud1 = ucell_c.atoms[T1].taud[I1];

		for(int j=0; j<cmd_neigh.nlist[i]; j++)
		{
            int T2 = ucell_c.iat2it[cmd_neigh.list[i][j]];
            int I2 = ucell_c.iat2ia[cmd_neigh.list[i][j]];
            Vector3<double> taud2 = ucell_c.atoms[T2].taud[I2];

			Vector3<double> tempd = cmd_neigh.cell_periodic(taud1,taud2);
            Vector3<double> temp;
            Mathzone::Direct_to_Cartesian(
			tempd.x, tempd.y, tempd.z,
			ucell_c.latvec.e11, ucell_c.latvec.e12, ucell_c.latvec.e13,
			ucell_c.latvec.e21, ucell_c.latvec.e22, ucell_c.latvec.e23,
			ucell_c.latvec.e31, ucell_c.latvec.e32, ucell_c.latvec.e33,
			temp.x, temp.y, temp.z);

            double distance = temp.norm()*ucell_c.lat0;

			if(distance <= INPUT.mdp.rcut_lj)
			{
				potential += LJ_energy(distance) - LJ_energy(INPUT.mdp.rcut_lj);
				force[i] = force[i] + LJ_force(distance, temp*ucell_c.lat0);
			}
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