#include "numerical_vna_lm.h"

Numerical_Vna_Lm::Numerical_Vna_Lm()
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

	r_radial = new double[1];
	k_radial = new double[1];
	rab = new double[1];
	psi = new double[1];
	psir = new double[1];
	psik = new double[1];
	
}

Numerical_Vna_Lm::~Numerical_Vna_Lm()
{
	delete[] r_radial;
	delete[] k_radial;
	delete[] rab;
	delete[] psi;
	delete[] psir;
	delete[] psik;
}

void Numerical_Vna_Lm::set_vna_proj_info
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
	const double &lat0
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
	
	this->cal_kradial();
	this->norm_test();
	this->plot();

	return;
}

//use Sbt_new
void Numerical_Vna_Lm::cal_kradial(void)
{
	assert( this->nr > 0);
	double *jl = new double[nr];
	double *integrated_func = new double[nr];

	const double pref = sqrt( 2.0 / PI );
	//Sbt method
	
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
void Numerical_Vna_Lm::norm_test(void)const
{
//	TITLE(ofs_onscaling, "Numerical_Vna_Lm", "norm_test");
	//double asum_r = 0.0;
	//double asum_k = 0.0;

	// note here psir = psi * r
	double *f = new double[nr];
	for(int ir=0; ir<nr; ir++)
	{
		f[ir] = this->psir[ir] * this->psir[ir];
	}

	double sumr = 0.0;
	//double sumk = 0.0;

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

void Numerical_Vna_Lm::plot(void)const
{
	TITLE("Numerical_Vna_Lm","plot");
	
	string orbital_type;
	switch( this->angular_momentum_l )
	{
		case 0: orbital_type = "s"; break;
		case 1: orbital_type = "p"; break;
		case 2: orbital_type = "d"; break;
		case 3: orbital_type = "f"; break;
		case 4: orbital_type = "g"; break;
		default: WARNING_QUIT("Numerical_Vna_Lm::plot","Please check in functoin.");
	}

	if(MY_RANK==0)
	{
		stringstream ssr, ssk, ssru;
		ssr << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << index_chi+1 << "-vna-proj-r.dat";

		ssk << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << index_chi+1 << "-vna-proj-k.dat";

		ofstream ofsr(ssr.str().c_str());
		ofstream ofsk(ssk.str().c_str());

		if (!ofsk || !ofsr)
		{
			WARNING("Numerical_Vna_Lm : plot", "Can't open files !");
		}

		for (int i = 0; i < this->nr; i++)
		{
			ofsr << this->r_radial[i] << " " << psi[i] << endl;
		}

		for (int i = 0; i < this->nk; i++)
		{
			ofsk << this->k_radial[i] << " " << psik[i] << endl;
		}

		ofsr.close();
		ofsk.close();
	}
	return;
}
