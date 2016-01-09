//=========================================================
//AUTHOR : liaochen, mohan
//DATE : 2008-11-12
//=========================================================
#include "numerical_orbital_lm.h"

Numerical_Orbital_Lm::Numerical_Orbital_Lm()
{
	label = "";
	index_atom_type = 0;
	angular_momentum_l = 0;
	index_chi = 0;

	nr=1;
	nk=1;

	rcut=0.0;
	kcut=0.0;
	dk=0.0;

	nr_uniform = 1;
	dr_uniform = -1.0;
	zty = 0.0;

	psi_uniform = new double[1];
	dpsi_uniform = new double[1];
	r_radial = new double[1];
	k_radial = new double[1];
	rab = new double[1];
	psi = new double[1];
	psir = new double[1];
	psik = new double[1];
	
}

Numerical_Orbital_Lm::~Numerical_Orbital_Lm()
{
	if(test_deconstructor)
	{
		cout << " ~Numerical_Orbital_Lm()" << endl;
	}
	delete[] psi_uniform;
	delete[] dpsi_uniform;
	delete[] r_radial;
	delete[] k_radial;
	delete[] rab;
	delete[] psi;
	delete[] psir;
	delete[] psik;
}

void Numerical_Orbital_Lm::set_orbital_info
(
 	const string &label_in,
	const int &index_atom_type_in,
	const int &angular_momentum_l_in,
	const int &index_chi_in,
	const int &nr_in,
	const double *rab_in,
	const double *r_radial_in,
	const double *psi_in,
	const int &nk_in,
	const double &dk_in,
	const double &lat0,
	const double &dr_uniform_in
)
{
    this->label = label_in;
    this->index_atom_type = index_atom_type_in;
    this->angular_momentum_l = angular_momentum_l_in;
    this->index_chi = index_chi_in;
    assert(nr_in>=2);
    assert(nr_in<10000);
    assert(nr%2!=0);
    this->nr = nr_in;
    assert(r_radial_in[nr-1]>0.0);
    assert(r_radial_in[nr-1]<50);
    assert(lat0 > 0);
    this->rcut = r_radial_in[nr-1];
    assert(nk_in>1);
    assert(nk_in<10000);
    this->nk = nk_in;
    assert(nk%2!=0);
    assert(dk_in>0);
    this->dk = dk_in;
	this->dr_uniform=dr_uniform_in;
	
	delete[] r_radial;
	delete[] k_radial;
	delete[] rab;
	delete[] psi;
	delete[] psir;
	delete[] psik;

	r_radial = new double[nr];
	k_radial = new double[nk];
	rab = new double[nr];
	psi = new double[nr];
	psir = new double[nr];
	psik = new double[nk];

	/***********************************************************
	be careful! LiaoChen modify on 2010/4/21
	************************************************************/
//	this->dk = PI / rcut / 2.0;
//	this->nk = this->nr;

	for (int ir = 0; ir < nr; ir++)
	{
		this->r_radial[ir] = r_radial_in[ir];
		this->rab[ir] = rab_in[ir];
		this->psi[ir] = psi_in[ir];
		this->psir[ir] = psi[ir] * r_radial[ir]; //mohan 2010-04-19
	}

	for (int ik = 0; ik < nk; ik++)
	{
		this->k_radial[ik] = ik * this->dk;
	}
	this->kcut = (nk-1) * this->dk;

	//liaochen modify on 2010/4/7
	//we do SBT on regular mesh
	//so we first generate psi_uniform first
	//we put uniform in ahead of cal_kradial
	
	bool uni = true;
	if (uni)
	{
		this->extra_uniform(dr_uniform);
	}
	else
	{
		this->use_uniform(dr_uniform);
	}

	this->cal_kradial();
	this->norm_test();
	this->plot();

	return;
}

void Numerical_Orbital_Lm::extra_uniform(const double &dr_uniform_in)
{
	timer::tick("NOrbital_Lm", "extra_uniform");
	
	//---------------------------------------------
	// set the dr, fixed by liaochen.
	// calculate the number of radial mesh points.
	//---------------------------------------------
	assert(dr_uniform>0.0);
	this->dr_uniform = dr_uniform_in;
	this->nr_uniform = static_cast<int>(rcut/dr_uniform) + 10;
	
	delete[] this->psi_uniform;
	this->psi_uniform = new double[nr_uniform];
	ZEROS (this->psi_uniform, nr_uniform);

	// do interpolation here to make grid more dense
	for (int ir = 0; ir < this->nr_uniform; ir++)
	{
		this->psi_uniform[ir] = Mathzone_Add1::Uni_RadialF(this->psi, this->nr, this->rab[0], ir * dr_uniform); 
//    	this->psi_uniform[ir] = Mathzone::Polynomial_Interpolation(this->psi, this->nr, this->rab[0], ir * dr_uniform); 
    }
	
	//----------------------------------------------	 
	// calculate the dpsi_uniform
	//----------------------------------------------	 
	delete [] this->dpsi_uniform;	
	this->dpsi_uniform = new double[this->nr_uniform];

	double* y2 = new double[nr];

	//--------------------------------------------------------------------------
	// old code to calculate the derivate dpsi/dr, 
	// has problem that the derivatives of orbitals oscillate a lot
	// around r=0
	//--------------------------------------------------------------------------
	//Mathzone_Add1::SplineD2 (r_radial, psi, nr, 100000.0, 100000.0, y2);
	//double yp1=(this->psi[1]-this->psi[0])/this->r_radial[1];
	//cout<<"psi0="<<"  "<<this->psi[0]<<"  "<<"psi1="<<"  "<<this->psi[1]<<"  "<<"r1="<<"  "<<this->r_radial[1]<<endl; 
	//cout<<"yp1="<<"  "<<yp1<<endl;
	//Mathzone_Add1::SplineD2 (r_radial, psi, nr, yp1, 0.0, y2);
	

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// new code developed by pengfei.
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	switch( this->angular_momentum_l ) // added by pengfei 13-8-8 different l has different  boundary conditions 
	{
		case 0: Mathzone_Add1::SplineD2 (r_radial, psi, nr, 0.0, 0.0, y2);; break;
		case 1: Mathzone_Add1::SplineD2 (r_radial, psi, nr, 100000.0, 100000.0, y2);; break;
		case 2: Mathzone_Add1::SplineD2 (r_radial, psi, nr, 0.0, 0.0, y2);; break;
		case 3: Mathzone_Add1::SplineD2 (r_radial, psi, nr, 0.0, 0.0, y2);; break;
		case 4: Mathzone_Add1::SplineD2 (r_radial, psi, nr, 0.0, 0.0, y2);; break;
		default: 
		  ofs_running << " The angular momentum should not be larger than 4 (g orbitals)." << endl;
		  exit(0);
	}

	//Mathzone_Add1::SplineD2 (r_radial, psi, nr, 0.0, 0.0, y2);
	//cout<<"angular_momentum_l="<<"  "<<this->angular_momentum_l<<endl;
	//for (int i=0; i<nr; i++)
	//{
	//     cout<<r_radial[i]<<"  "<<y2[i]<<endl;
	//}
	//Method 1
	//	Mathzone_Add1::Uni_Deriv_Phi (psi_uniform, nr_uniform, dr_uniform, 1, dpsi_uniform);
	//	Mathzone_Add1::Uni_Deriv_Phi (psi_uniform, nr_uniform, dr_uniform, 2, ddpsi_uniform);




	double* rad = new double[nr_uniform];
	for (int ir = 0; ir < nr_uniform; ir++)
	{
		rad[ir] = ir*dr_uniform;
	}

	//	Mathzone_Add1::SplineD2 (rad, psi_uniform, nr_uniform, 0.0, 0.0, ddpsi_uniform);
	double* tmp = new double[nr_uniform];
	Mathzone_Add1::Cubic_Spline_Interpolation(r_radial, psi, y2, 
			nr, rad, nr_uniform, tmp, dpsi_uniform );

	// calculate zty
	// liaochen add 2010-08
	Mathzone_Add1::Uni_Deriv_Phi (this->psi_uniform, this->nr_uniform, dr_uniform, angular_momentum_l, tmp);
	this->zty = tmp[0]/Mathzone_Add1::factorial (angular_momentum_l);

	delete [] y2;
	delete [] rad;
	delete [] tmp;
	timer::tick("NOrbital_Lm", "extra_uniform");
}

void Numerical_Orbital_Lm::use_uniform(const double &dr_uniform_in)
{
	assert(dr_uniform_in>0.0);
	this->dr_uniform = dr_uniform_in;
	// for save: +10, because in real space interpolation,
	// there may be "one grid point" more than the cutoff.
	this->nr_uniform = static_cast<int>(rcut/dr_uniform)+10;

	delete[] this->psi_uniform;
	this->psi_uniform = new double[nr_uniform];
	ZEROS(psi_uniform, nr_uniform);

	string orbital_type;
	switch( this->angular_momentum_l )
	{
		case 0: orbital_type = "s"; break;
		case 1: orbital_type = "p"; break;
		case 2: orbital_type = "d"; break;
		case 3: orbital_type = "f"; break;
		case 4: orbital_type = "g"; break;
	}
		
	cout << "===========================================================" << endl;
	for(int i=0; i<nr_uniform; i++)
	{
		this->psi_uniform[i] = 
			Mathzone_Add1::Uni_RadialF(this->psi, this->nr, this->rab[0], i*dr_uniform); 
	}
	
	delete [] this->dpsi_uniform;
	this->dpsi_uniform = new double[nr_uniform];
	Mathzone_Add1::Uni_Deriv_Phi (psi_uniform, nr_uniform, dr_uniform, 1, dpsi_uniform);
	
	if(MY_RANK==0)
	{
		stringstream ss;
		ss << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << ".ORBITAL_NOR_uniform.txt";

		ofstream ofs(ss.str().c_str());

		for(int i=0; i<nr_uniform; i++)
		{
			ofs << setw(15) << i*dr_uniform << setw(20) << psi_uniform[i] << endl;
		}
		ofs.close();
	}

	return;
}

//liaochen modify on 2010/4/7
//use Sbt_new
void Numerical_Orbital_Lm::cal_kradial(void)
{
	assert( this->nr > 0);
	assert( this->nr_uniform > 0);
	double *jl = new double[nr];
	double *integrated_func = new double[nr];

	const double pref = sqrt( 2.0 / PI );
	//Sbt method
	
	/*
	double* rad = new double[nr_uniform];
	for (int ir = 0; ir < nr_uniform; ir++) 
	{
		rad[ir] = dr_uniform * ir;
	}
	
	//liaochen add
	Mathzone_Add1::Sbt_new (3, angular_momentum_l, 
							k_radial, dk, nk, 
							rad, dr_uniform, nr_uniform, 
							psi_uniform, 0, this->psik);
	
	for (int ik = 0; ik < nk; ik++) this->psik[ik] *= (pref*k_radial[ik]);
	delete [] rad;
	*/
	
	//integration directly
	for (int ik = 0; ik < nk; ik++)
	{
		Mathzone::Spherical_Bessel(
				this->nr, 
				this->r_radial, 
				this->k_radial[ik], 
				this->angular_momentum_l, 
				jl);

		for (int ir = 0; ir < nr; ir++)
		{
			integrated_func[ir] = this->psir[ir] * this->r_radial[ir] * jl[ir];
		}

		Mathzone::Simpson_Integral(
				this->nr, 
				integrated_func, 
				this->rab, 
				this->psik[ik]);

		this->psik[ik] *= ( pref * k_radial[ik]);
	}
		
	delete[] integrated_func;
	delete[] jl;
}

//===============================================
//FOUND LOCAL VARIABLE
//asum : integral of psi*psi in whole space
//===============================================
void Numerical_Orbital_Lm::norm_test(void)const
{
//	TITLE(ofs_onscaling, "Numerical_Orbital_Lm", "norm_test");
	double asum_r = 0.0;
	double asum_k = 0.0;

	// note here psir = psi * r
	double *f = new double[nr];
	for(int ir=0; ir<nr; ir++)
	{
		f[ir] = this->psir[ir] * this->psir[ir];
	}

	double sumr = 0.0;
	double sumk = 0.0;

	Mathzone::Simpson_Integral(this->nr, f, this->rab, sumr);

	delete[] f;
	f = new double[nk];
	for(int ik=0; ik<nk; ik++)
	{
		f[ik] = this->psik[ik] * this->psik[ik];
	}

//	Mathzone::Simpson_Integral(this->nk, f, this->k_radial, sumk);
	
	//means nothing.
	//ofs_running << setw(12) << sumk << endl;

	delete[] f;
	return;
}

void Numerical_Orbital_Lm::plot(void)const
{
	TITLE("Numerical_Orbital_Lm","plot");
	
	string orbital_type;
	switch( this->angular_momentum_l )
	{
		case 0: orbital_type = "s"; break;
		case 1: orbital_type = "p"; break;
		case 2: orbital_type = "d"; break;
		case 3: orbital_type = "f"; break;
		case 4: orbital_type = "g"; break;
		default: WARNING_QUIT("Numerical_Orbital_Lm::plot","Please check in functoin.");
	}

	if(MY_RANK==0)
	{
		stringstream ssr, ssk, ssru ,ssdru; // 2013-08-10 pengfei
		ssr << global_out_dir << this->label << "/"
			<< this->label << "-"<< orbital_type << index_chi+1 << "-orbital-r.dat";

		ssk << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << index_chi+1 << "-orbital-k.dat";

		ssru << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << index_chi+1 << "-orbital-ru.dat";

		ssdru << global_out_dir << this->label << "/"  // 2013-08-10 pengfei
			<< this->label << "-" << orbital_type << index_chi+1 << "-orbital-dru.dat";

		ofstream ofsr(ssr.str().c_str());
		ofstream ofsk(ssk.str().c_str());
		ofstream ofsru(ssru.str().c_str());
		ofstream ofsdru(ssdru.str().c_str()); // 2013-08-10 pengfei

		if (!ofsk || !ofsr || !ofsru || !ofsdru) // 2013-08-10 pengfei
		{
			WARNING("Numerical_Orbital_Lm : plot", "Can't open files !");
		}

		for (int i = 0; i < this->nr; i++)
		{
			ofsr << this->r_radial[i] << " " << psi[i] << endl;
		}

		for (int i = 0; i < this->nk; i++)
		{
			ofsk << this->k_radial[i] << " " << psik[i] << endl;
		}

		for (int i = 0; i < this->nr_uniform; i++)
		{
			ofsru << this->dr_uniform * i << " " << psi_uniform[i] << endl;
		}

		for (int i = 0; i < this->nr_uniform; i++)
		{
			ofsdru << this->dr_uniform * i << " " << dpsi_uniform[i] << endl;// output dphi/dr 2013-08-10  pengfei
		}


		ofsr.close();
		ofsk.close();
		ofsru.close();
		ofsdru.close(); // 13-08-10 pengfei
	}
	return;
}
