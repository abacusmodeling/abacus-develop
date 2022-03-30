#include "pseudo_nc.h"

pseudo_nc::pseudo_nc()
{
// pseudo_h
	els = new std::string[1];
	lchi = new int[1];
	oc = new double[1];
	jjj = new double[1];
	jchi = new double[1];
	nn = new int[1];
	has_so = false;
	zv = 0;

// pseudo local parts
	vloc_at =  new double[1];

// pseudo_atom
	r = new double[1];
	rab = new double[1];
	rho_at = new double[1];
	rho_atc = new double[1];

// pseudo_nc
	lll = new int[1];
}

pseudo_nc::~pseudo_nc()
{
// pseudo_h
	delete[] els;
	delete[] lchi;
	delete[] oc;
	delete[] jjj;
	delete[] jchi;
	delete[] nn;

// pseudo local parts
	delete[] vloc_at;

// pseudo_atom
	delete[] r;
	delete[] rab;
	delete[] rho_at;
	delete[] rho_atc;

// pseudo_nc
	delete[] lll;
}



//---------------------------------------------------------------------
void pseudo_nc::set_pseudo_nc(const Pseudopot_upf &upf)
{
	ModuleBase::TITLE("pseudo_nc","set_pseudo_nc");

	// call subroutines
	this->set_pseudo_h(upf);
	this->set_pseudo_atom(upf);
	this->set_pseudo_vl(upf);


	delete[] lll;
	lll = new int[nbeta];
	assert(lll != 0);

	for (int i = 0;i < nbeta;i++)
	{
		lll[i] = upf.lll[i];
	}

	kkbeta = 0;

	for (int nb = 0;nb < nbeta;nb++)
	{
		kkbeta = (upf.kkbeta[nb] > kkbeta) ? upf.kkbeta[nb] : kkbeta;
	}

	betar.create(upf.beta.nr, upf.beta.nc);
	
//	OUT("betar.nr",upf.beta.nr); // nbeta
//	OUT("betar.nc",upf.beta.nc); // mesh

	dion.create(nbeta, nbeta);

	betar = upf.beta;

	dion = upf.dion;

	nh = 0;

	for (int nb = 0; nb < nbeta;nb++)
	{
		nh += 2 * lll [nb] + 1;
	} 

	return;
} // end subroutine set_pseudo_upf


void pseudo_nc::print_pseudo_nc(std::ofstream &ofs)
{
	print_pseudo_vl(ofs);
	ofs << "\n pseudo_nc : ";
	ofs << "\n kkbeta	" << kkbeta;
	ofs << "\n nh  " << nh;
	output::printr1_d(ofs, " lll : ", lll, nbeta);
	output::printrm(ofs, " betar : ", betar);
	output::printrm(ofs, " dion : ", dion);
	ofs << "\n ----------------------";
}


void pseudo_nc::set_pseudo_h(const Pseudopot_upf &upf)
{
	ModuleBase::TITLE("pseudo_nc::set_pseudo_h");
	// set pseudopotential for each atom type
	// by using the Unified Pseudopotential Format

	this->nv = upf.nv;// UPF file version number
	this->psd = upf.psd;
	this->pp_type = upf.pp_type;

	this->tvanp = upf.tvanp;// if USPP
	this->nlcc = upf.nlcc;// Non linear core corrections( bool ?)
	
	this->xc_func = upf.xc_func;

	this->zv = upf.zp;
	this->etotps = upf.etotps;
	this->ecutwfc = upf.ecutwfc;
	this->ecutrho = upf.ecutrho;

	this->lmax = upf.lmax;
	this->mesh = upf.mesh;
	
	// mohan update 2021-02-22
	//  max number of points in the atomic radial mesh
	int ndmx = 2000; 
	if (this->mesh > ndmx)
	{
		std::cout << "\n set_pseudo_h, too many grid points,";
	}

	this->nchi = upf.nwfc;
	this->nbeta = upf.nbeta;

	delete[] els;
	this->els = new std::string[nchi];
	assert(this->els != 0);

	delete[] lchi;
	this->lchi = new int[this->nchi];
	assert(this->lchi != 0);

	delete[] oc;
	this->oc = new double[nchi];
	assert(this->oc != 0);

	for (int i = 0;i < nchi;i++)
	{
		this->els[i] = upf.els[i];
		this->lchi[i] = upf.lchi[i];
		this->oc[i] = upf.oc[i];
	}

	delete[] jjj;
	this->jjj = new double[nbeta];
	assert(this->jjj != 0);

	delete[] nn;
	this->nn = new int[nchi];
	assert(this->nn != 0);

	delete[] jchi;
	this->jchi = new double[nchi];
	assert(this->jchi != 0);

	this->has_so = upf.has_so;//added by zhengdy-soc

	if (this->has_so)
	{ 
		for(int i=0;i < nchi;i++)
		{
			this->nn[i] = upf.nn[i];
			this->jchi[i] = upf.jchi[i];
		}
		for(int i=0;i < upf.nbeta;i++)
		{
			this->jjj[i]  = upf.jjj [i];
		}
	}
	else
	{
		for (int i=0; i<nchi; i++)
		{
			this->nn[i] = 0;
			this->jchi[i] = 0;
		}
		for (int i=0; i<upf.nbeta; i++)
		{
			this->jjj[i]  = 0;
		}
	} 

	/*
		jchi = new double[nwfc];
		if (upf.has_so){
			for(i=0;i<upf.nwfc;i++){
				jchi[i] = upf.jchi[i];
			}
		}else{
			for(i=0;i<upf.nwfc;i++){
				jchi[i] = 0;
			}
		} // endif
	*/

	return;
} // end subroutine set_pseudo_upf


void pseudo_nc::set_pseudo_atom(const Pseudopot_upf &upf)
{
	ModuleBase::TITLE("pseudo_nc","set_pseudo_atom");

	// mohan 2009-12-15
	// mohan update again 2011-05-23, 
	// in order to calculate more accurate Vna.
	this->rcut = GlobalV::PSEUDORCUT;//(a.u.);
	
// remember to update here if you need it.
//	rcut = 25.0; 

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"PAO radial cut off (Bohr)",rcut);
	if(rcut <= 0.0)
	{
		ModuleBase::WARNING_QUIT("pseudo_atom::set_pseudo_atom","PAO rcut<=0.0");
	}

	chi.create(nchi, mesh);

	delete[] r;
	r = new double[mesh];
	assert(r != 0);
	ModuleBase::GlobalFunc::ZEROS(r, mesh);

	delete[] rab;
	rab = new double[mesh];
	assert(rab != 0);
	ModuleBase::GlobalFunc::ZEROS(rab, mesh);

	delete[] rho_at;
	rho_at  = new double[mesh];
	assert(rho_at != 0);
	ModuleBase::GlobalFunc::ZEROS(rho_at,mesh);

	delete[] rho_atc;
	rho_atc = new double[mesh];
	assert(rho_atc != 0);
	ModuleBase::GlobalFunc::ZEROS(rho_atc, mesh);

	for (int i = 0;i < nchi;i++)
	{
		for (int j = 0; j < mesh; j++)
		{
			chi(i, j) = upf.chi(i, j);
		}
	}

	for (int i = 0;i < mesh;i++)
	{
		r  [i]     = upf.r  [i];
		rab[i]     = upf.rab[i];
		rho_at [i] = upf.rho_at [i];
	}

	// nlcc: non-linear core corrections
	// rho_atc: core atomic charge
	if (nlcc)
	{
		for (int i = 0;i < mesh;i++)
		{
			rho_atc[i] = upf.rho_atc[i];
		}
	}
	else
	{
		for (int i = 0;i < upf.mesh;i++)
		{
			rho_atc[i] = 0.0;
		}
	} // end if


	bool br = false;

	this->msh = 0;

	for (int ir = 0;ir < mesh;ir++)
	{
		if (r [ir] > rcut)
		{
			msh = ir + 1;
			br = true;
			break;
		}
	}

	if (br)
	{
		// force msh to be odd for simpson integration
		msh = 2 * (int)((msh + 1) / 2) - 1;	// 5
	}
	else
	{
		msh = mesh ;
	}

	return;
} // end subroutine set_pseudo



void pseudo_nc::set_pseudo_vl(const Pseudopot_upf &upf)
{
	ModuleBase::TITLE("pseudo_nc","set_pseudo_vl");

	assert(mesh>0);//mohan add 2021-05-01

	delete[] vloc_at;
	vloc_at = new double[mesh];
	assert(vloc_at != 0);

	for (int i = 0;i < mesh;i++)
	{
		vloc_at[i] = upf.vloc[i];
	}

	return;
} 


void pseudo_nc::print_pseudo_atom(std::ofstream &ofs)
{
	print_pseudo_h(ofs);
	ofs << "\n pseudo_atom : ";
	ofs << "\n msh	" << msh;
//	ofs	<< "\n nchi	" << nchi;
	output::printr1_d(ofs, " r : ", r, mesh);
	output::printr1_d(ofs, " rab : ", rab, mesh);
	output::printr1_d(ofs, " rho_atc : ", rho_atc, mesh);
	output::printr1_d(ofs, " rho_at : ", rho_at, mesh);
	output::printr1_d(ofs," jchi : ", jchi, nchi);
	output::printrm(ofs, " chi : ", chi);
	ofs << "\n ----------------------";
}


void pseudo_nc::print_pseudo_vl(std::ofstream &ofs)
{
	ofs << "\n pseudo_vl:";
	print_pseudo_atom(ofs);
	output::printr1_d(ofs, "vloc_at : ", vloc_at, mesh);
	ofs << "\n ----------------------------------- ";
}

void pseudo_nc::print_pseudo_h(std::ofstream &ofs)
{
    ofs << "\n pseudo_info :";
    ofs << "\n nv       " << nv;
    ofs << "\n psd  " << psd;
    ofs << "\n pp_type  " << pp_type;
    ofs << "\n tvanp    " << tvanp;
    ofs << "\n nlcc " << nlcc;
    ofs << "\n dft  " << xc_func;
    ofs << "\n zv       " << zv;
    ofs << "\n etotps   " << etotps;
    ofs << "\n ecutwfc " << ecutwfc;
    ofs << "\n ecutrho " << ecutrho;
    ofs << "\n lmax " << lmax;
    ofs << "\n mesh " << mesh;
    ofs << "\n nchi " << nchi;
    ofs << "\n nbeta    " << nbeta;
//  out.printr1_d(ofs," els: ", els, nchi);
    output::printr1_d(ofs, " lchi: ", lchi, nchi);
    output::printr1_d(ofs, " oc: ", oc, nchi);
    ofs << "\n ----------------------";
}

