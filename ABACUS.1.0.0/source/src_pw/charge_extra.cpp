#include "charge_extra.h"
#include "../src_pw/tools.h"
#include "global.h"
#include "../src_lcao/global_fp.h"

Charge_Extra::Charge_Extra()
{
	init_rho = false;
	//xiaohui add 2014-05-07, for first-order extrapolation
	this->delta_rho1 = new double*[NSPIN];
	this->delta_rho2 = new double*[NSPIN];
	this->delta_rho = new double*[NSPIN];
	//xiaohui add 2014-05-10, for second-order extrapolation
	this->delta_rho3 = new double*[NSPIN];

	for(int is=0; is<NSPIN; is++)
	{
		delta_rho1[is] = new double[pw.nrxx];
		delta_rho2[is] = new double[pw.nrxx];
		delta_rho[is] = new double[pw.nrxx];
		//xiaohui add 2014-05-10, for second-order extrapolation
		delta_rho3[is] = new double[pw.nrxx];

		ZEROS(delta_rho1[is], pw.nrxx);
		ZEROS(delta_rho2[is], pw.nrxx);
		ZEROS(delta_rho[is], pw.nrxx);
		//xiaohui add 2014-05-10, for second-order extrapolation
		ZEROS(delta_rho3[is], pw.nrxx);
	}
	pos_old1 = new double[1];
	pos_old2 = new double[1];
	pos_now = new double[1];
	pos_next = new double[1];
	alpha = 0.0;
	beta = 0.0;
}

Charge_Extra::~Charge_Extra()
{
	if(init_rho)
	{
		for(int i=0; i<dim; i++)
		{
			for(int is=0; is<NSPIN; is++)
			{
				delete[] rho_ion[i][is];
			}
			delete[] rho_ion[i];
		}	
	}
	//xiaohui add 2014-05-07, for first-order extrapolation
	for(int is=0; is<NSPIN; is++)
	{
		delete[] delta_rho1[is];
		delete[] delta_rho2[is];
		delete[] delta_rho[is];
		//xiaohui add 2014-05-10, for second-order extrapolation
		delete[] delta_rho3[is];
         }	
	 delete[] delta_rho1;
	 delete[] delta_rho2;
 	 delete[] delta_rho;
	 //xiaohui add 2014-05-10, for second-order extrapolation
	 delete[] delta_rho3;

	//xiaohui add 2014-05-10, for second-order extrapolation
	delete[] pos_old1;
	delete[] pos_old2;
	delete[] pos_now;
	delete[] pos_next;
}

void Charge_Extra::allocate(void)
{
	TITLE("Charge_Extra","allocate");

	// 1: first order extrapolation.
	// 2: second order extrapolation.
	//xiaohui modify 2014-05-11, for second-order extrapolation	
	//if(pot.extra_pot != 1 && pot.extra_pot != 2)
	//{
	//	return;
	//}

	this->dim = 0;//xiaohui modify 2015-02-01
	//xiaohui add 2014-05-10, for second-order extrapolation	
	pos_dim = ucell.nat * 3;

	delete[] this->pos_old1;
	delete[] this->pos_old2;
	delete[] this->pos_now;
	delete[] this->pos_next;

	this->pos_old1 = new double[pos_dim];
	this->pos_old2 = new double[pos_dim];
	this->pos_now = new double[pos_dim];
	this->pos_next = new double[pos_dim];

	ZEROS(pos_old1, pos_dim);
	ZEROS(pos_old2, pos_dim);
	ZEROS(pos_now, pos_dim);
	ZEROS(pos_next, pos_dim);

	/* xiaohui modify 2014-05-12, test
	//xiaohui add 2014-05-11, test
	delete[] this->pos_old1;
	this->pos_old1 = new double[pos_dim];
	ZEROS(pos_old1, pos_dim);
	
	delete[] this->pos_old2;
	this->pos_old2 = new double[pos_dim];
	ZEROS(pos_old2, pos_dim);

	//delete[] this->pos_now;
	//this->pos_now = new double[pos_dim];
	//ZEROS(pos_now, pos_dim);

	//delete[] this->pos_next;
	//this->pos_next = new double[pos_dim];
	//ZEROS(pos_next, pos_dim);
	*/

	//cout<<"here pos_dim is: "<<pos_dim<<endl;
	
	//cout<<"It is ok here. 2014-05-11"<<endl;

	if(init_rho)
	{
		WARNING_QUIT("Charge_Extra::allocate","rho_ion has been allocated, pls check.");
	}

	rho_ion = new double**[dim];

	for(int i=0; i<dim; i++)
	{
		rho_ion[i] = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			rho_ion[i][is] = new double[pw.nrxx];
			// first value from charge density.
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				rho_ion[i][is][ir] = chr.rho[is][ir];	
			}
		}
	}	

	init_rho = true;

	Memory::record("charge_extra","rho_ion",dim*NSPIN*pw.nrxx,"double");

	return;
}


void Charge_Extra::record_rho()
{
	TITLE("Charge_Extra","record_rho");

}

void Charge_Extra::extrapolate_charge()
{
    TITLE("Charge_Extra","extrapolate_charge");
    // charge density expolation:
    // pot_order = 0 copy the old potential(nothing is done);
    // pot_order = 3 substrate old atomic charge density and
    // sum the new.
    // If the dynamics is done, this routine extrapolates also
    // the difference between the scf and the atomic one.
    // pot_order = 1 first order extrapolation:
    // rho(t+dt) = 2*rho(t) - rho(t-dt)
    // pot_order = 2 second order extrapolation:
    // rho(t+dt) = rho+
    // alpha*(rho(t) - rho(t-dt))
    // + beta*( rho(t-dt) - rho(t-2*dt) )
	// just use atomic charge.
/*xiaohui modify 2015-02-01
	if(pot.extra_pot == 0)
	{
		// do nothing
		//xiaohui add 2014-05-10, setup structure factor for extra_pot=0
		ofs_running << " Setup the structure factor in plane wave basis." << endl;
		pw.setup_structure_factor();
	}
	else if(pot.extra_pot == 1)
	{
		// 1st order
		for(int is=0; is<NSPIN; is++)
		{
			double *rho_tmp = new double[pw.nrxx];
			ZEROS(rho_tmp, pw.nrxx);
				
			// (1) save the current charge density.
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				rho_tmp[ir] = chr.rho[is][ir];	
			} 
			
			// (2) prepare charge density for the next ion iteration.
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				chr.rho[is][ir] = 2.0*rho_tmp[ir] - this->rho_ion[0][is][ir];
			}

			// (3) save the current charge density for next step.
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				this->rho_ion[0][is][ir] = rho_tmp[ir];
			}

			delete[] rho_tmp;
		}
	}
	else if(pot.extra_pot == 2)
	{
		// 2 nd order		
	}
	else if(pot.extra_pot == 3)
	{
		// start from atomic charge density.
		// worst method ever!
		chr.atomic_rho(NSPIN, chr.rho);
	}
xiaohui modify 2015-02-01*/
	//else if(pot.extra_pot == 4)
	if(pot.extra_pot == "dm")//xiaohui modify 2015-02-01
	{
		//if(LOCAL_BASIS!=4 || LINEAR_SCALING!=1) xiaohui modify 2013-09-01
		if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-01
		{
			WARNING_QUIT("Charge_Extra","extrapolate_charge");
		}
		else
		{
			pw.setup_structure_factor();
			// should not do this after grid_technique is done!.
//			for(int is=0; is<NSPIN; is++)
//			{
//				ZEROS(chr.rho[is], pw.nrxx);
//			}
// 			UHM.GG.cal_rho();
		}
	}
	//xiaohui add 2014-05-03, "atomic" extrapolation
	else if(pot.extra_pot == "atomic")
	{
		//xiaohui add 2014-05-03, "atomic" extrapolation
		//cout<<"atomic extrapolation begin"<<endl;

		double** rho_atom_old = new double*[NSPIN];
		double** rho_atom_new = new double*[NSPIN];
		//double** delta_rho = new double*[NSPIN];

		for(int is=0; is<NSPIN; is++)
		{
			rho_atom_old[is] = new double[pw.nrxx];
			rho_atom_new[is] = new double[pw.nrxx];
			//delta_rho[is] = new double[pw.nrxx];

			ZEROS(rho_atom_old[is], pw.nrxx);
			ZEROS(rho_atom_new[is], pw.nrxx);
			//ZEROS(delta_rho[is], pw.nrxx);
		}
		chr.atomic_rho(NSPIN,rho_atom_old);
		for(int is=0; is<NSPIN; is++)
		{
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				delta_rho[is][ir] = chr.rho[is][ir] - rho_atom_old[is][ir];
			}
		}

		//xiaohui add 'OUT_LEVEL', 2015-09-16
		if(OUT_LEVEL != "m") ofs_running << " Setup the structure factor in plane wave basis." << endl;
		pw.setup_structure_factor();

		//xiaohui add 2014-05-03
		chr.atomic_rho(NSPIN,rho_atom_new);
		for(int is=0; is<NSPIN; is++)
                {
                	for(int ir=0; ir<pw.nrxx; ir++)
                        {
                        	chr.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
                        }
                }
		for(int is=0; is<NSPIN; is++)
               	{
                       	delete[] rho_atom_old[is];
			delete[] rho_atom_new[is];
			//delete[] delta_rho[is];
                }	
               	delete[] rho_atom_old;
		delete[] rho_atom_new;
		//delete[] delta_rho;
		//cout<<"atomic extrapolation finish"<<endl;
		
	}
	//xiaohui add 2014-05-07, "first-order" extrapolation
	else if(pot.extra_pot == "first-order")
	{
		//xiaohui add 2014-05-07, "first-order" extrapolation
		//cout<<"first-order extrapolation begin"<<endl;

		double** rho_atom_old = new double*[NSPIN];
		double** rho_atom_new = new double*[NSPIN];
		//this->delta_rho1 = new double*[NSPIN];
		//this->delta_rho2 = new double*[NSPIN];
		//this->delta_rho = new double*[NSPIN];

		for(int is=0; is<NSPIN; is++)
		{
			rho_atom_old[is] = new double[pw.nrxx];
			rho_atom_new[is] = new double[pw.nrxx];
			//delta_rho1[is] = new double[pw.nrxx];
			//delta_rho2[is] = new double[pw.nrxx];
			//delta_rho[is] = new double[pw.nrxx];

			ZEROS(rho_atom_old[is], pw.nrxx);
			ZEROS(rho_atom_new[is], pw.nrxx);
			//ZEROS(delta_rho1[is], pw.nrxx);
			//ZEROS(delta_rho2[is], pw.nrxx);
			//ZEROS(delta_rho[is], pw.nrxx);
		}
		chr.atomic_rho(NSPIN,rho_atom_old);
		for(int is=0; is<NSPIN; is++)
		{
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				delta_rho2[is][ir] = delta_rho1[is][ir];
				delta_rho1[is][ir] = chr.rho[is][ir] - rho_atom_old[is][ir];
				delta_rho[is][ir] = 2*delta_rho1[is][ir] - delta_rho2[is][ir];
			}
		}

		//xiaohui add "OUT_LEVEL", 2015-09-16
		if(OUT_LEVEL != "m") ofs_running << " Setup the structure factor in plane wave basis." << endl;
		pw.setup_structure_factor();

		//xiaohui add 2014-05-07
		chr.atomic_rho(NSPIN,rho_atom_new);
		for(int is=0; is<NSPIN; is++)
                {
                	for(int ir=0; ir<pw.nrxx; ir++)
                        {
				if(istep == 1)
				{
                        		chr.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
				}
				else
				{
                        		chr.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
				}
                        }
                }
		for(int is=0; is<NSPIN; is++)
               	{
                       	delete[] rho_atom_old[is];
			delete[] rho_atom_new[is];
                }	
               	delete[] rho_atom_old;
		delete[] rho_atom_new;
		//cout<<"first-order extrapolation finish"<<endl;
		
	}

	//xiaohui add 2014-05-10, "second-order" extrapolation
	else if(pot.extra_pot == "second-order")
	{
		//xiaohui add 2014-05-10, "second-order" extrapolation
		//cout<<"second-order extrapolation begin"<<endl;

		double** rho_atom_old = new double*[NSPIN];
		double** rho_atom_new = new double*[NSPIN];
		//this->delta_rho1 = new double*[NSPIN];
		//this->delta_rho2 = new double*[NSPIN];
		//this->delta_rho = new double*[NSPIN];

		for(int is=0; is<NSPIN; is++)
		{
			rho_atom_old[is] = new double[pw.nrxx];
			rho_atom_new[is] = new double[pw.nrxx];
			//delta_rho1[is] = new double[pw.nrxx];
			//delta_rho2[is] = new double[pw.nrxx];
			//delta_rho[is] = new double[pw.nrxx];

			ZEROS(rho_atom_old[is], pw.nrxx);
			ZEROS(rho_atom_new[is], pw.nrxx);
			//ZEROS(delta_rho1[is], pw.nrxx);
			//ZEROS(delta_rho2[is], pw.nrxx);
			//ZEROS(delta_rho[is], pw.nrxx);
		}

		chr.atomic_rho(NSPIN,rho_atom_old);

		find_alpha_and_beta();
		//cout<<"alpha= "<<alpha<<endl;
		//cout<<"beta= "<<beta<<endl;

		for(int is=0; is<NSPIN; is++)
		{
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				delta_rho3[is][ir] = delta_rho2[is][ir];
				delta_rho2[is][ir] = delta_rho1[is][ir];
				delta_rho1[is][ir] = chr.rho[is][ir] - rho_atom_old[is][ir];
				delta_rho[is][ir] = delta_rho1[is][ir] + 
							alpha * (delta_rho1[is][ir] - delta_rho2[is][ir]) +
								beta * (delta_rho2[is][ir] - delta_rho3[is][ir]);
			}
		}

		//xiaohui add 'OUT_LEVEL', 2015-09-16
		if(OUT_LEVEL != "m") ofs_running << " Setup the structure factor in plane wave basis." << endl;
		pw.setup_structure_factor();

		//xiaohui add 2014-05-07
		chr.atomic_rho(NSPIN,rho_atom_new);
		for(int is=0; is<NSPIN; is++)
                {
                	for(int ir=0; ir<pw.nrxx; ir++)
                        {
				if(istep == 1)
				{
                        		chr.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
				}
				else if(istep == 2)
				{
					delta_rho[is][ir] = 2*delta_rho1[is][ir] - delta_rho2[is][ir];
					chr.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
				}
				else
				{
                        		chr.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
				}
                        }
                }
		for(int is=0; is<NSPIN; is++)
               	{
                       	delete[] rho_atom_old[is];
			delete[] rho_atom_new[is];
                }	
               	delete[] rho_atom_old;
		delete[] rho_atom_new;
		//cout<<"second-order extrapolation finish"<<endl;
	}
	else
	{
		WARNING_QUIT("potential::init_pot","extra_pot parameter is wrong!");
	}


    return;
}

void Charge_Extra::find_alpha_and_beta()
{
	double a11,a12,a21,a22;
	double b1,b2;
	double detA;

	a11 = 0.0;
	a12 = 0.0;
	a21 = 0.0;
	a22 = 0.0;
	b1 = 0.0;
	b2 = 0.0;
	detA = 0.0;

	if(istep >= 3)
	{
		int iat=0;
		for(int it = 0;it < ucell.ntype;it++)
		{
			//Atom* atom = &ucell.atoms[it];
			for(int ia =0;ia< ucell.atoms[it].na;ia++)
			{
				a11 += (pos_now[3*iat  ] - pos_old1[3*iat  ]) * (pos_now[3*iat  ] - pos_old1[3*iat  ]) + 
					(pos_now[3*iat+1] - pos_old1[3*iat+1]) * (pos_now[3*iat+1] - pos_old1[3*iat+1]) + 
					(pos_now[3*iat+2] - pos_old1[3*iat+2]) * (pos_now[3*iat+2] - pos_old1[3*iat+2]);

				a12 += (pos_now[3*iat  ] - pos_old1[3*iat  ]) * (pos_old1[3*iat  ] - pos_old2[3*iat  ]) + 
					(pos_now[3*iat+1] - pos_old1[3*iat+1]) * (pos_old1[3*iat+1] - pos_old2[3*iat+1]) + 
					(pos_now[3*iat+2] - pos_old1[3*iat+2]) * (pos_old1[3*iat+2] - pos_old2[3*iat+2]);

				a22 += (pos_old1[3*iat  ] - pos_old2[3*iat  ]) * (pos_old1[3*iat  ] - pos_old2[3*iat  ]) + 
					(pos_old1[3*iat+1] - pos_old2[3*iat+1]) * (pos_old1[3*iat+1] - pos_old2[3*iat+1]) + 
					(pos_old1[3*iat+2] - pos_old2[3*iat+2]) * (pos_old1[3*iat+2] - pos_old2[3*iat+2]);

				b1 += -((pos_now[3*iat  ] - pos_next[3*iat  ]) * (pos_now[3*iat  ] - pos_old1[3*iat  ]) +
					(pos_now[3*iat+1] - pos_next[3*iat+1]) * (pos_now[3*iat+1] - pos_old1[3*iat+1]) +
					(pos_now[3*iat+2] - pos_next[3*iat+2]) * (pos_now[3*iat+2] - pos_old1[3*iat+2]));

				b2 += -((pos_now[3*iat  ] - pos_next[3*iat  ]) * (pos_old1[3*iat  ] - pos_old2[3*iat  ]) +
					(pos_now[3*iat+1] - pos_next[3*iat+1]) * (pos_old1[3*iat+1] - pos_old2[3*iat+1]) +
					(pos_now[3*iat+2] - pos_next[3*iat+2]) * (pos_old1[3*iat+2] - pos_old2[3*iat+2]));

				iat++;
			}
		}
	}
	a21 = a12;
	detA = a11 * a22 - a21 * a12;

	//xiaohui add 2015-05-12, test
	//cout<<"a11= "<<a11<<endl;
	//cout<<"a12= "<<a12<<endl;
	//cout<<"a21= "<<a21<<endl;
	//cout<<"a22= "<<a22<<endl;
	//cout<<"b1= "<<b1<<endl;
	//cout<<"b2= "<<b2<<endl;
	//cout<<"detA= "<<detA<<endl;

	if(detA > 1.0E-12)
	{
		alpha = (b1 * a22 - b2 * a12) / detA;
		beta = (b2 * a11 - b1 * a21) / detA;
		if(abs(alpha) >10) alpha=1.0;
		if(abs(beta)>10) beta=0.0;
	}
	else
	{
		//xiaohui modify 2014-08-04
		//alpha = 1.0;
		alpha = 0.0;
		beta = 0.0;
		if(a11 != 0)
		{
			alpha = b1 / a11 ;
		}
		if(abs(alpha) >10) alpha=1.0;
	}
	return;
}
