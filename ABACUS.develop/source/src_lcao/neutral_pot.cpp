#include "neutral_pot.h"
#include "../src_pw/tools.h"
#include "../src_pw/global.h"
#include "../src_global/poission.h"

Neutral_Pot::Neutral_Pot()
{
	vna = new double[1];
	vna_u = new double[1];

	// for Vna projector
	nchi = new int[1];

	lmax = -1;
	mesh_u = 0;
}

Neutral_Pot::~Neutral_Pot()
{
	delete[] vna;
	delete[] vna_u;

	// for Vna projectors
	delete[] nchi;
}

void Neutral_Pot::setup_Vna(const int &it)
{
	TITLE("Neutral_Pot","setup_Vna");

	ofs_running << "\n SETUP THE 1D NEUTRAL POTENTIAL VNA FOR " << ucell.atoms[it].label << endl;
	//cout << " setup the neutral potential Vna." << endl;

	//==================================================
	// setup the local part of pseudopotentials.
	//==================================================
	Atom* atom = &ucell.atoms[it];
	const int msh = atom->msh;
	double* vloc = atom->vloc_at;
	double* r = atom->r;
	stringstream ss;
	ss << global_out_dir << atom->label << "/" << atom->label << "_vl_pseudo_r.dat";
	this->output_1dv(ss.str(),msh,vloc,r);

	//==================================================
	// setup the one dimensional atomic charge density.
	//==================================================
	double* rhoatm = new double[msh];
	double* vatom = new double[msh];
	ZEROS(rhoatm, msh);
	ZEROS(vatom, msh);

	for(int ir=0; ir<msh; ir++)
	{
		double r2 = r[ir] * r[ir];
		rhoatm[ir] = atom->rho_at[ir] / FOUR_PI / r2;
	}

	double charge = 0.0;
	Mathzone::Simpson_Integral(atom->msh,atom->rho_at,atom->rab,charge);

	OUT(ofs_running,"charge from pseudo file",charge);

	//---------------------------------------------
	// need to scale the charge if 
	// the charge read from the pseudopotential
	// is not equal to the actul valence electron
	// number
	//---------------------------------------------
	assert(charge!=0.0);
	double scale=1.0;
	if(charge!=atom->zv)
	{
		OUT(ofs_running,"charge should be",atom->zv);
		scale = atom->zv/charge;
	}

	for(int ir=0; ir<msh; ++ir)
	{
		rhoatm[ir] *= scale;
	}


	Poission::SolPoissonEq(rhoatm,r,msh,vatom);
	delete[] rhoatm;

	stringstream ssc;
	ssc << global_out_dir << atom->label << "/" << atom->label << "_v_atom.dat";
	this->output_1dv(ssc.str(),msh,vatom,r);

	//==================================
	// calculate the neutral potential.
	//==================================
	delete[] vna;
	vna = new double[msh];
	ZEROS(vna, msh);
	
	for(int ir=0; ir<msh; ir++)
	{
		vna[ir] = vatom[ir] + vloc[ir];
		//vna[ir] = vloc[ir]; //tmp by mohan
	}
	delete[] vatom;

	stringstream ssn;
	ssn << global_out_dir << atom->label << "/" << atom->label << "_vna.dat";
	this->output_1dv(ssn.str(),msh,vna,r);

	return;
}

void Neutral_Pot::output_1dv(const string &fn, const int &msh, const double* target, const double* r)const
{
	if(MY_RANK!=0) return;
	ofstream ofs(fn.c_str());
	for(int ir=0; ir<msh; ir++)
	{
		ofs << r[ir] << " " << target[ir] << endl;
	}
	ofs.close();
	return;
}

void Neutral_Pot::uniform_Vna(const int &it, const double &dr_in )
{
	assert( vna!=0 );
	assert( dr_in > 0.0 );

	this->dr = dr_in;
	Atom* atom = &ucell.atoms[it];
	this->rcut = atom->r[ atom->msh-1 ];

	double threshold = 1.0e-5; //tmp by mohan
	//if threshold goes down to 1.0e-5, oscillation may happen.
	for(int ir=0; ir< atom->msh; ir++)
	{
		if( abs( vna[ir] ) < threshold )
		{
			this->rcut = atom->r[ir];
			break;
		}
	}

	// mohan note: 2011-05-23,
	// if rcut > 15, check pseudo_atom.cpp, 
	// the readin radius must be larger.
	assert(rcut <= 15);

	OUT(ofs_running,"radius cutoff of Vna (Bohr)",rcut);
	OUT(ofs_running,"threshold of Vna",threshold);

	this->nr = (int)(rcut / dr_in)+5;//mohan fix bug 2011-06-03

	delete[] vna_u;
	this->vna_u = new double[nr];
	ZEROS(vna_u, nr);

	double* r = new double[nr];
	ZEROS(r, nr);

    for (int ir = 0; ir < this->nr; ir++)
    {
        r[ir] = ir * this->dr;
		// mohan fix serious bug 2012-01-11,
		// the dimension of atom->r should be the 3rd parameter.
        //this->vna_u[ir] = Mathzone::Polynomial_Interpolation_xy( atom->r, vna, this->nr, r[ir]);
        this->vna_u[ir] = Mathzone::Polynomial_Interpolation_xy( atom->r, vna, atom->msh, r[ir]);
    }

	stringstream ssu;
	ssu << global_out_dir << atom->label << "/" << atom->label << "_vna_uniform.dat";
	this->output_1dv(ssu.str(),nr,vna_u,r);

	delete[] r;

	return;
}



void Neutral_Pot::init_proj(Numerical_Orbital &Phi, const double &dr_uniform)
{
	cout << " ==> init_proj of Vna";
	// set the lmax we used
	this->lmax = Phi.getLmax();

	// set the total number of Chi, which equals the Vna projector now.
	this->total_nchi;
	for(int L=0; L<=lmax; ++L)
	{
		this->total_nchi += Phi.getNchi(L);
	}

	cout << " lmax=" << lmax << " total_nchi=" << total_nchi << endl;
	
	delete[] nchi;
	this->nchi = new int[lmax+1];	
	ZEROS( nchi, lmax+1);

	for(int L=0; L<=lmax; ++L)
	{
		nchi[L] = Phi.getNchi(L);
	}

	// initialize the Chi = Gram-Schemit ( VnaPsi )
	cout << " nr for Vna is " << this->nr << endl;
	cout << " rcut fo Vna is " << this->rcut << endl;
	cout << " Rcut for Orbital is " << Phi.getRcut() << endl;


	double radius = std::min( this->rcut, Phi.getRcut() );
	assert(radius>0.0);
	assert(dr_uniform>0.0);
	this->mesh_u = static_cast<int>(radius / dr_uniform);
	if(mesh_u%2==0) ++mesh_u;
	cout << " mesh_u = " << mesh_u << endl;
	
	// calculate the radius function Phi * Vna
	double** VnaPhi = new double*[total_nchi];
	for(int ic=0; ic<total_nchi; ++ic)
	{
		VnaPhi[ic] = new double[mesh_u]; 
	}


	// set Vna*Phi and then Gram-Schmit them...
	this->set_vnaphi(VnaPhi, Phi, dr_uniform);
	
	
	
	for(int ic=0; ic<total_nchi; ++ic)
	{
		delete[] VnaPhi[ic];
	}
	delete[] VnaPhi;
	return;
}







void Neutral_Pot::set_vnaphi(double** VnaPhi, Numerical_Orbital &Phi, const double &dr_uniform)
{

	return;
/*
	//multiply Vna by atomic orbitals.
	int ic_index=0;
	for(int L=0; L<= this->lmax; ++L)
	{
		for(int N=0; N< this->nchi[L]; ++N)
		{
			for(int ir=0; ir< this->mesh_u; ++ir)
			{
				VnaPhi[ic_index][ir] = this->vna_u[ir] * Phi.PhiLN(L,N).psi_uniform[ir];
			}
			++ic_index;
		}
	}

	// 
	for(int L=0; L<= this->lmax; ++L)
	{
		for(int N=0; N< this->nchi[L]; ++N)
		{
			double overlap = 0.0;
			for(int ir=0; ir< this->mesh_u; ++ir)
			{
	//			orbr[i] = VnaPhi[ic_index][ir] * ((double)ir*ORB.dr_uniform); 
			}
			++ic_index;	
		}
	}


	//-------------------------------------
	// Schmit-Gramit orthogonalization.
	//-------------------------------------
	double **new_orb = new double*[this->total_nchi];
	for(int i=0; i<total_nchi; ++i)
	{
		ZEROS(new_orb[i], mesh_u);
	}

	double* rab = new double[mesh_u];
	for(int ir=0; ir<mesh_u; ++ir)
	{
		rab[ir] = dr_uniform;
	}
	
	//--------------------------------------------------
	// Implement the detail of Schmit-Gramit methods.
	//--------------------------------------------------
	ic_index = 0;
	for(int L=0; L<= this->lmax; ++L)
	{
		for(int N1=0; N1< this->nchi[L]; ++N1)
		{
			if(N1==0)
			{
				for(int ir=0; ir< this->mesh_u; ++ir)
				{
					new_orb[ic_index][ir] = vna_u[ir] * Phi.PhiLN(L,N1).psi_uniform[ir]; 
				}
			}

			//------------------------------------------------
			// calculate the overlap < Phi_new | Vna | Phi >
			// and < Phi_new | Vna | Phi_new >
			//------------------------------------------------
			else if(N1>0)
			{
				const int NN = nchi[L] - 1;
				// must orthogonal to the previous multi-zeta orbitals.
				for(int N_pre=0; N_pre<NN; ++N_pre)
				{
					double overlap = 0.0;
					double* f = new double[mesh_u];

					// Vna * Phi * r^2  
					for(int ir=0; ir<mesh_u; ++ir)
					{
						f[ir] = VnaPhi[ic_index][ir] * pow( ir*ORB.dr_uniform, 2.0 ) ;
					}
					Mathzone::Simpson_Integral(this->mesh_u, f, rab, overlap);

					// ddot < Phi_new | Vna | Phi >
					for(int ir=0; ir<mesh_u; ++ir)
					{

					}	


					delete[] f;
				}
			}
			++ic_index;
		}
	}


	delete[] rab;

	// delete the overlap.
	for(int i=0; i<total_nchi; ++i)
	{
		delete[] new_orb[i];
	}
	delete[] new_orb;

*/
	return;
}
