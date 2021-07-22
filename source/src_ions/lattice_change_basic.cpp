#include "lattice_change_basic.h"
#include "../src_pw/global.h"

int Lattice_Change_Basic::dim=0;
bool Lattice_Change_Basic::converged=true;
double Lattice_Change_Basic::largest_grad=0.0;
int Lattice_Change_Basic::update_iter=0;
int Lattice_Change_Basic::istep=0;

double Lattice_Change_Basic::ediff=0.0;
double Lattice_Change_Basic::etot=0.0;
double Lattice_Change_Basic::etot_p=0.0;

//double Lattice_Change_Basic::lattice_change_ini = 0.5; // default is 0.5
double Lattice_Change_Basic::lattice_change_ini = 0.01; // default is 0.5

int Lattice_Change_Basic::out_stru=0;

void Lattice_Change_Basic::setup_gradient(double* lat, double *grad, matrix &stress)
{
	TITLE("Lattice_Change_Basic","setup_gradient");
	
	if(INPUT.fixed_axes == "volume")
	{
		double stress_aver = (stress(0,0) + stress(1,1) + stress(2,2))/3.0;
		stress(0,0) = stress(0,0) - stress_aver;
		stress(1,1) = stress(1,1) - stress_aver;
		stress(2,2) = stress(2,2) - stress_aver;
	}
	
	lat[0] = ucell.latvec.e11 * ucell.lat0;   lat[1] = ucell.latvec.e12 * ucell.lat0; lat[2] = ucell.latvec.e13 * ucell.lat0; 
	lat[3] = ucell.latvec.e21 * ucell.lat0;   lat[4] = ucell.latvec.e22 * ucell.lat0; lat[5] = ucell.latvec.e23 * ucell.lat0;
	lat[6] = ucell.latvec.e31 * ucell.lat0;   lat[7] = ucell.latvec.e32 * ucell.lat0; lat[8] = ucell.latvec.e33 * ucell.lat0;
	
	if(ucell.lc[0] == 1)
	{
		grad[0] = -(lat[0]*stress(0,0) + lat[1]*stress(1,0) + lat[2]* stress(2,0));
		grad[1] = -(lat[0]*stress(0,1) + lat[1]*stress(1,1) + lat[2]* stress(2,1));
		grad[2] = -(lat[0]*stress(0,2) + lat[1]*stress(1,2) + lat[2]* stress(2,2));
	}
	if(ucell.lc[1] == 1)
	{
		grad[3] = -(lat[3]*stress(0,0) + lat[4]*stress(1,0) + lat[5]* stress(2,0));
		grad[4] = -(lat[3]*stress(0,1) + lat[4]*stress(1,1) + lat[5]* stress(2,1));
		grad[5] = -(lat[3]*stress(0,2) + lat[4]*stress(1,2) + lat[5]* stress(2,2));
	}
	if(ucell.lc[2] == 1)
	{
		grad[6] = -(lat[6]*stress(0,0) + lat[7]*stress(1,0) + lat[8]* stress(2,0));
		grad[7] = -(lat[6]*stress(0,1) + lat[7]*stress(1,1) + lat[8]* stress(2,1));
		grad[8] = -(lat[6]*stress(0,2) + lat[7]*stress(1,2) + lat[8]* stress(2,2));
	}

        //grad[0] = -stress(0,0);   grad[1] = -stress(0,1);  grad[2] = -stress(0,2);
        //grad[3] = -stress(1,0);   grad[4] = -stress(1,1);  grad[5] = -stress(1,2);
        //grad[6] = -stress(2,0);   grad[7] = -stress(2,1);  grad[8] = -stress(2,2);
	
	return;
}

void Lattice_Change_Basic::change_lattice(double *move, double *lat)
{
	TITLE("Lattice_Change_Basic","change_lattice");

	assert(move!=NULL);
	assert(lat!=NULL);

/*
	cout<<" LATTICE CONSTANT  OLD:"<<endl;
	cout<<" "<<setprecision(12)<<ucell.latvec.e11<<"   "<<ucell.latvec.e12<<"   "<<ucell.latvec.e13<<endl;
	cout<<" "<<setprecision(12)<<ucell.latvec.e21<<"   "<<ucell.latvec.e22<<"   "<<ucell.latvec.e23<<endl;
	cout<<" "<<setprecision(12)<<ucell.latvec.e31<<"   "<<ucell.latvec.e32<<"   "<<ucell.latvec.e33<<endl;
*/
	
	if(ucell.lc[0] != 0)
	{
		ucell.latvec.e11 = (move[0] + lat[0])/ucell.lat0;
		ucell.latvec.e12 = (move[1] + lat[1])/ucell.lat0;
		ucell.latvec.e13 = (move[2] + lat[2])/ucell.lat0;
	}
	if(ucell.lc[1] !=0)
	{
		ucell.latvec.e21 = (move[3] + lat[3])/ucell.lat0;
		ucell.latvec.e22 = (move[4] + lat[4])/ucell.lat0;
		ucell.latvec.e23 = (move[5] + lat[5])/ucell.lat0;
	}
	if(ucell.lc[2] !=0)
	{
		ucell.latvec.e31 = (move[6] + lat[6])/ucell.lat0;
		ucell.latvec.e32 = (move[7] + lat[7])/ucell.lat0;
		ucell.latvec.e33 = (move[8] + lat[8])/ucell.lat0;
	}
	
	ucell.a1.x = ucell.latvec.e11;
	ucell.a1.y = ucell.latvec.e12;
	ucell.a1.z = ucell.latvec.e13;
	ucell.a2.x = ucell.latvec.e21;
	ucell.a2.y = ucell.latvec.e22;
	ucell.a2.z = ucell.latvec.e23;
	ucell.a3.x = ucell.latvec.e31;
	ucell.a3.y = ucell.latvec.e32;
	ucell.a3.z = ucell.latvec.e33;

	ucell.omega = abs( ucell.latvec.Det() ) * ucell.lat0 * ucell.lat0 * ucell.lat0;
	
	ucell.GT = ucell.latvec.Inverse();
	ucell.G  = ucell.GT.Transpose();
	ucell.GGT = ucell.G * ucell.GT;
	ucell.invGGT = ucell.GGT.Inverse();

#ifdef __MPI
    // distribute lattice vectors.
    Parallel_Common::bcast_double(ucell.latvec.e11 );
    Parallel_Common::bcast_double(ucell.latvec.e12 );
    Parallel_Common::bcast_double(ucell.latvec.e13 );
    Parallel_Common::bcast_double(ucell.latvec.e21 );
    Parallel_Common::bcast_double(ucell.latvec.e22 );
    Parallel_Common::bcast_double(ucell.latvec.e23 );
    Parallel_Common::bcast_double(ucell.latvec.e31 );
    Parallel_Common::bcast_double(ucell.latvec.e32 );
    Parallel_Common::bcast_double(ucell.latvec.e33 );
	
    // distribute lattice vectors.
    Parallel_Common::bcast_double( ucell.a1.x );
    Parallel_Common::bcast_double( ucell.a1.y );
    Parallel_Common::bcast_double( ucell.a1.z );
    Parallel_Common::bcast_double( ucell.a2.x );
    Parallel_Common::bcast_double( ucell.a2.y );
    Parallel_Common::bcast_double( ucell.a2.z );
    Parallel_Common::bcast_double( ucell.a3.x );
    Parallel_Common::bcast_double( ucell.a3.y );
    Parallel_Common::bcast_double( ucell.a3.z );
#endif
/*
        cout<<" LATTICE CONSTANT NEW: "<<endl;
        cout<<" "<<setprecision(12)<<ucell.latvec.e11<<"   "<<ucell.latvec.e12<<"   "<<ucell.latvec.e13<<endl;
        cout<<" "<<setprecision(12)<<ucell.latvec.e21<<"   "<<ucell.latvec.e22<<"   "<<ucell.latvec.e23<<endl;
        cout<<" "<<setprecision(12)<<ucell.latvec.e31<<"   "<<ucell.latvec.e32<<"   "<<ucell.latvec.e33<<endl;
*/

	
	return;
}

void Lattice_Change_Basic::check_converged(matrix &stress, double *grad)
{
	TITLE("Lattice_Change_Basic","check_converged");

	Lattice_Change_Basic::largest_grad = 0.0;
	double stress_ii_max = 0.0;
	
	if(ucell.lc[0] == 1 && ucell.lc[1] == 1 && ucell.lc[2] == 1)
	{
		for(int i=0;i<3;i++)
		{
			if(stress_ii_max < abs(stress(i,i))) stress_ii_max = abs(stress(i,i));
			for(int j=0;j<3;j++)
			{
				if(Lattice_Change_Basic::largest_grad < abs(stress(i,j)))
				{
					Lattice_Change_Basic::largest_grad = abs(stress(i,j));
				}
			}
		}
	}
	else
	{
		for(int i=0; i<9; i++)
		{
			if(Lattice_Change_Basic::largest_grad < abs(grad[i]))
			{
				Lattice_Change_Basic::largest_grad = abs(grad[i]);
			}
		}
	}

	double unit_transform = 0.0;
	unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
	Lattice_Change_Basic::largest_grad = Lattice_Change_Basic::largest_grad * unit_transform;
	stress_ii_max = stress_ii_max * unit_transform;
	
	if(Lattice_Change_Basic::largest_grad == 0.0)
	{
		ofs_running << " largest stress is 0, no movement is possible." << endl;
		ofs_running << " it may converged, otherwise no movement of lattice parameters is allowed." << endl;
		Lattice_Change_Basic::converged = true;
	}
	else if(ucell.lc[0] == 1 && ucell.lc[1] == 1 && ucell.lc[2] == 1)
	{
		//if(Lattice_Change_Basic::largest_grad < STRESS_THR)
		if(Lattice_Change_Basic::largest_grad < STRESS_THR && stress_ii_max < STRESS_THR)
		{
			ofs_running << "\n Lattice relaxation is converged!" << endl;
			ofs_running << "\n Largest gradient is = " << largest_grad << endl;
			Lattice_Change_Basic::converged = true;
			++ Lattice_Change_Basic::update_iter;
		}
		else
		{
			 ofs_running << "\n Lattice relaxation is not converged yet (threshold is "<< STRESS_THR << ")" << endl;
			 Lattice_Change_Basic::converged = false;
		}
	}
	else
	{
		/*for(int i=0; i<9; i++)
		{
			cout<<"i= "<<i<<" "<<grad[i]<<endl;
		}*/
		if(Lattice_Change_Basic::largest_grad < 10 * STRESS_THR)
		{
			ofs_running << "\n Lattice relaxation is converged!" << endl;
			ofs_running << "\n Largest gradient is = " << largest_grad << endl;
			Lattice_Change_Basic::converged = true;
			++ Lattice_Change_Basic::update_iter;
		}
		else
		{
			 ofs_running << "\n Lattice relaxation is not converged yet (threshold is "<< STRESS_THR << ")" << endl;
			 Lattice_Change_Basic::converged = false;
		}
	}

	return;		
}


void Lattice_Change_Basic::terminate(void)
{
	TITLE("Lattice_Change_Basic","terminate");
	if(Lattice_Change_Basic::converged)
	{
		ofs_running << " end of lattice optimization"<<endl;
		OUT(ofs_running,"istep", Lattice_Change_Basic::istep);
		OUT(ofs_running,"update iteration", Lattice_Change_Basic::update_iter);
		/*
		ofs_running<<"Saving the approximate inverse hessian"<<endl;
		ofstream hess("hess.out");
		for(int i=0;i<dim;i++)
		{
			for(int j=0;j<dim;j++)
			{
				hess << inv_hess(i,j);
			}
		}
		hess.close();
		*/
	}
	else
	{
		ofs_running<<" the maximum number of steps has been reached." << endl;
		ofs_running<<" end of lattice optimization." << endl;
	}

	return;
}

void Lattice_Change_Basic::setup_etot(const double &energy_in, const bool judgement)
{
	if(Lattice_Change_Basic::istep==1)
	{
		// p == previous
		Lattice_Change_Basic::etot_p = energy_in;
		Lattice_Change_Basic::etot   = energy_in;
		ediff = etot - etot_p;
	}
	else
	{
		if(judgement) 
		{
			Lattice_Change_Basic::etot = energy_in;
			if(Lattice_Change_Basic::etot_p > etot)
			{
				ediff = etot - etot_p;
				Lattice_Change_Basic::etot_p = etot;
			}
			else
			{
				// this step will not be accepted
				ediff = 0.0;
			}
		}	
		else // for bfgs
		{
			Lattice_Change_Basic::etot_p = etot;
			Lattice_Change_Basic::etot = energy_in;
			ediff = etot - etot_p;
		}
	}

	return;	
}







