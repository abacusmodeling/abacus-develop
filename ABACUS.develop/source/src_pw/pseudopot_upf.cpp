#include "../src_io/output.h"
#include "pseudopot_upf.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring> // Peize Lin fix bug about strcpy 2016-08-02


using namespace std;

Pseudopot_upf::Pseudopot_upf()
{
	this->els = new string[1];
	this->lchi = new int[1];
	this->oc = new double[1];

	this->r = new double[1];
	this->rab = new double[1];
	this->vloc = new double[1];

	this->kkbeta = new int[1];
	this->lll = new int[1];

	this->rho_at = new double[1];
	this->rho_atc = new double[1];

	this->nn = new int[1];//zhengdy-soc
	this->jchi = new double[1];
	this->jjj = new double[1];

	functional_error = 0;//xiaohui add 2015-03-24
}

Pseudopot_upf::~Pseudopot_upf()
{
	delete [] els; 
	delete [] lchi;
	delete [] oc;

	delete [] r;    //mesh_1
	delete [] rab;  //mesh_2
	delete [] vloc;  //local_1

	delete [] kkbeta; // nl_1
	delete [] lll; // nl_2

	delete [] rho_at;// psrhoatom_1
	delete [] rho_atc;

	delete [] nn;
	delete [] jjj;
	delete [] jchi;
}

int Pseudopot_upf::init_pseudo_reader(const string &fn)
{
    TITLE("Pseudopot_upf","init");
    // First check if this pseudo-potential has spin-orbit information
    ifstream ifs(fn.c_str(), ios::in);

	// can't find the file.
	if (!ifs)
    {
        return 1;
    }

    if(global_pseudo_type=="auto") //zws
	{
		set_pseudo_type(fn);
	}

	// read in the .UPF type of pseudopotentials
	if(global_pseudo_type=="upf")
	{
		int info = read_pseudo_upf(ifs);
		return info;
	}
	// read in the .vwr type of pseudopotentials
	else if(global_pseudo_type=="vwr")
	{
		int info = read_pseudo_vwr(ifs);
		return info;
	}
	else if(global_pseudo_type=="upf201")
	{
		int info = read_pseudo_upf201(ifs);
		return info;
	}

	return 0;
}


//----------------------------------------------------------
// setting the type of the pseudopotential file
//----------------------------------------------------------
int Pseudopot_upf::set_pseudo_type(const string &fn) //zws add
{
    ifstream pptype_ifs(fn.c_str(), ios::in);
    string dummy;
	string strversion;

	if (pptype_ifs.good())
	{
		getline(pptype_ifs,dummy);

		stringstream wdsstream(dummy);
		getline(wdsstream,strversion,'"');
		getline(wdsstream,strversion,'"');

		if ( trim(strversion) == "2.0.1" )
		{
			global_pseudo_type = "upf201";
		}
		else
		{
			global_pseudo_type = "upf";
		}
	}
	return 0;
}

string& Pseudopot_upf::trim(string &in_str)
{
    static const string deltri = " \t" ; // delete tab or space
    string::size_type position = in_str.find_first_of(deltri, position);
    if (position == string::npos)
	{
        return in_str;
	}
    return trim(in_str.erase(position, 1) );
}

string Pseudopot_upf::trimend(string &in_str)
{
    const string &deltri =" \t" ;
    string::size_type position = in_str.find_last_not_of(deltri)+1;
    string tmpstr=in_str.erase(position);
    return tmpstr.erase(0,tmpstr.find_first_not_of(deltri));
} //zws



//  read pseudopot_upf potential "upf" in the Unified
//  Pseudopot_upfpotential Format
int Pseudopot_upf::read_pseudo_upf(ifstream &ifs)
{
	string dummy;
	int ierr = 0;
	this->has_so = false;

	//addinfo_loop
	ifs.rdstate();

	while (ifs.good())
	{
		ifs >> dummy;
		if(dummy=="<PP_ADDINFO>")
		{
			ierr = 1;
			break;
		}
		ifs.rdstate();
	}

	if (ierr == 1) 
	{
		has_so = true;
	}
	
	// Search for Header
	// This version doesn't use the new routine SCAN_BEGIN
	// because this search must set extra flags for
	// compatibility with other pp format reading

	ierr = 0;

	ifs.clear();
	ifs.seekg(0);
	ifs.rdstate();

	// header_loop:
	while (ifs.good())
	{
		ifs >> dummy;
		if(dummy=="<PP_HEADER>")
		{
			ierr = 1;
			read_pseudo_header(ifs);
			SCAN_END(ifs, "</PP_HEADER>");
			break;
		}
	}

	if (ierr == 0) 
	{
		// 2: something in pseudopotential file not match.
		return 2;
	}

	// Search for mesh information
	if ( SCAN_BEGIN(ifs, "<PP_MESH>") )
	{
		read_pseudo_mesh(ifs);
		SCAN_END(ifs, "</PP_MESH>");
	}

	// If  present, search for nlcc
	if (this->nlcc)
	{
		SCAN_BEGIN(ifs, "<PP_NLCC>"); 
		read_pseudo_nlcc(ifs);
		SCAN_END(ifs, "</PP_NLCC>");
	}

	// Search for Local potential
	SCAN_BEGIN(ifs, "<PP_LOCAL>");
	read_pseudo_local(ifs);
	SCAN_END(ifs, "</PP_LOCAL>");

	// Search for Nonlocal potential
	SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	read_pseudo_nl(ifs);
	SCAN_END(ifs, "</PP_NONLOCAL>");

	// Search for atomic wavefunctions
	SCAN_BEGIN(ifs, "<PP_PSWFC>");
	read_pseudo_pswfc(ifs);
	SCAN_END(ifs, "</PP_PSWFC>");

	// Search for atomic charge
	SCAN_BEGIN(ifs, "<PP_RHOATOM>");
	read_pseudo_rhoatom(ifs);
	SCAN_END(ifs, "</PP_RHOATOM>");

	// Search for add_info
	if (has_so)
	{
//		WARNING_QUIT("read_pseudo_upf","Can't read UPF containing spin-orbital term.");
		SCAN_BEGIN (ifs,"<PP_ADDINFO>");//added by zhengdy-soc
		read_pseudo_so (ifs);
		SCAN_END (ifs,"</PP_ADDINFO>");
	}

	ifs.clear();
	ifs.seekg(0);

	// return 0: read in sucessfully.
	return 0;
}// end subroutine read_pseudopot_upf


void Pseudopot_upf::read_pseudo_header(ifstream &ifs)
{
	READ_VALUE(ifs, this->nv);// Version number
	READ_VALUE(ifs, this->psd);// Element label

	// Type of pseudo : NC or US
	READ_VALUE(ifs, this->pp_type);
	if(pp_type=="US")
	{
		this->tvanp = true;
	}
	else if(pp_type=="NC")
	{
		this->tvanp = false;
	}
	else
	{
		// A bug here!!! can't quit together.
		cout << " pp_type=" << pp_type << endl;
		WARNING_QUIT("Pseudopot_upf::read_pseudo_header","unknown pseudo type");
	}

	// If use nlcc
	string nlc;
	READ_VALUE(ifs, nlc);

	if (nlc == "T")
	{
		this->nlcc = true;
	}
	else
	{
		this->nlcc = false;
	}

	// mohan modify 2009-12-15
	ifs >> dft[0] >> dft[1] >> dft[2];
	READ_VALUE(ifs, this->dft[3]);
	
	// dft functional enforced to modify
	// mohan add 2010-07-15
	if(DFT_FUNCTIONAL!="none")
	{
		/*xiaohui modify 2015-03-24
		dft[0] = DFT_FUNCTIONAL;
		dft[1] = DFT_FUNCTIONAL;
		dft[2] = DFT_FUNCTIONAL;
		dft[3] = DFT_FUNCTIONAL;
		xiaohui modify 2015-03-24*/

		//xiaohui add 2015-03-23
		string dft_functional;
		if(dft[1] == "PZ")
		{
			dft_functional = "lda";
		}
		else if(dft[1] == "PW")
		{
			dft_functional = "pbe";
		}

		if(dft_functional != DFT_FUNCTIONAL)
		{
			functional_error = 1;

			cout << " dft_functional readin is: " << DFT_FUNCTIONAL << endl;
			cout << " dft_functional in pseudopot file is: " << dft_functional << endl;
			ofs_warning << " dft_functional readin is: " << DFT_FUNCTIONAL << endl;
			ofs_warning << " dft_functional in pseudopot file is: " << dft_functional << endl;
			//WARNING_QUIT("Pseudopot_upf::read_pseudo_header","input xc functional does not match that in pseudopot file");
		}
	}
	
	READ_VALUE(ifs, this->zp);
	READ_VALUE(ifs, this->etotps);

	ifs >> this->ecutwfc >> this->ecutrho;
	ifs.ignore(75, '\n');

	READ_VALUE(ifs, this->lmax);
	READ_VALUE(ifs, this->mesh);

	ifs >> this->nwfc >> this->nbeta ;
	ifs.ignore(75, '\n');
	ifs.ignore(75, '\n');

	delete[] els;
	delete[] lchi;
	delete[] oc;
	this->els = new string[nwfc];
	this->lchi = new int[nwfc];
	this->oc = new double[nwfc];

	ZEROS(lchi, nwfc); // angular momentum of each orbital
	ZEROS(oc, nwfc);//occupation of each orbital

	for(int i=0;i<nwfc;i++)
	{
		ifs >> els[i] >> this->lchi[i] >> this->oc[i];
	}
	return;
}

void Pseudopot_upf::read_pseudo_mesh(ifstream &ifs)
{
	assert(mesh>0);

	delete[] r;
	delete[] rab;
	this->r = new double[mesh];
	this->rab = new double[mesh];
	ZEROS(r,mesh);
	ZEROS(rab,mesh);

	int ir = 0;

	if( SCAN_BEGIN(ifs, "<PP_R>", false) )
	{
		for (ir = 0;ir < mesh;ir++)
		{
			ifs >> this->r[ir];
		}
		SCAN_END(ifs, "</PP_R>");
	}

	if( SCAN_BEGIN(ifs, "<PP_RAB>", false) )
	{
		for (ir = 0;ir < mesh;ir++)
		{
			ifs >> this->rab[ir];
		}
		SCAN_END(ifs, "</PP_RAB>");
	}
	return;
}

void Pseudopot_upf::read_pseudo_nlcc(ifstream &ifs)
{
	assert(mesh>0);
	delete[] rho_atc;
	this->rho_atc = new double[mesh];
	ZEROS(rho_atc, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_atc[ir];
	}
	return;
}

void Pseudopot_upf::read_pseudo_local(ifstream &ifs)
{
	assert(mesh>0);
	delete[] vloc;
	this->vloc = new double[mesh];
	ZEROS(vloc, mesh);

	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->vloc[ir];
	}

	
	// for test
/*
	ofstream ofs("vloc_UPF.dat");
	for(int ir=0; ir<mesh; ++ir)
	{
		ofs << r[ir] << " " << vloc[ir] << endl;
	}
	ofs.close();
*/

	return;
}

void Pseudopot_upf::read_pseudo_nl(ifstream &ifs)
{
//	int nb, mb, n, ir, idum, ldum, lp, i, ikk;
	int nb, mb, ir, idum;

	if (nbeta == 0)
	{
		delete[] kkbeta;
		delete[] lll;
		this->kkbeta = new int[1];
		this->lll = new int[1];
		this->beta.create(1, 1);
		this->dion.create(1, 1);
		return;
	}
	else
	{
		delete[] kkbeta;
		delete[] lll;
		this->kkbeta = new int[nbeta]; 
		this->lll = new int[nbeta]; 
		this->beta.create(nbeta , mesh);
		this->dion.create(nbeta , nbeta);

		for(int i=0;i<nbeta;i++)
		{
			SCAN_BEGIN(ifs, "<PP_BETA>", false);
			ifs >> idum;
			READ_VALUE(ifs, this->lll[i]);// nl_1
			READ_VALUE(ifs, this->kkbeta[i]);// nl_2
			// number of mesh points for projectors

			for (ir=0;ir<kkbeta[i];ir++)
			{
				ifs >> this->beta(i, ir);// nl_3

				// --------- FOR TEST ---------
				// mohan test bug
				if(i>2)
				{
					beta(i,ir) = 0.0;
				}
				// --------- FOR TEST ---------
			}
			SCAN_END(ifs, "</PP_BETA>");
		}

		// DIJ
		SCAN_BEGIN(ifs, "<PP_DIJ>", false);
		READ_VALUE(ifs, this->nd);//nl_4
		for (int i=0; i<this->nd; i++)
		{
			double swap;
			ifs >> nb >> mb >> swap;
			nb--;
			mb--;
			assert( nb >= 0); // mohan add 2011-03-10
			assert( mb >= 0);
			this->dion(mb, nb) = swap;// nl_5
			this->dion(nb, mb) = swap;
		}
		SCAN_END(ifs, "</PP_DIJ>");

		// QIJ
		if (tvanp)
		{
			WARNING_QUIT("read_pseudo_nl","Ultra Soft Pseudopotential not available yet.");
		}
		else // not tvanp
		{
		}
	}
	return;
}

void Pseudopot_upf::read_pseudo_pswfc(ifstream &ifs)
{
	this->chi.create(this->nwfc, this->mesh);
	for (int i=0;i<nwfc;i++)
	{
		string OrbitalName;
		int BelongToL = 0;
		double occupation = 0.0;
		string dummy;
		ifs >> OrbitalName >> BelongToL >> occupation >> dummy;
		for (int ir = 0;ir < mesh;ir++)
		{
			ifs >> this->chi(i, ir);
		}
	}
	return;
}

void Pseudopot_upf::read_pseudo_rhoatom(ifstream &ifs)
{
	delete[] rho_at;
	this->rho_at = new double[mesh];
	ZEROS(rho_at, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_at[ir];
	}
	return;
}

void Pseudopot_upf::read_pseudo_so(ifstream &ifs)
{
       //read soc info from upf, added by zhengdy-soc
       if(!this->has_so) return;
       delete[] nn;
       delete[] jchi;
       delete[] jjj;
       this->nn = new int[nwfc];
       this->jchi = new double[nwfc];
       this->jjj = new double[nbeta];
       ZEROS(nn,nwfc);
       ZEROS(jchi,nwfc);
       ZEROS(jjj,nbeta);
       //RELWFC
       for(int nw=0;nw< nwfc;nw++)
       {
             ifs >> this->els[nw] >>this->nn[nw] >> this->lchi[nw] >> this->jchi[nw] >> this->oc[nw];
             if(this->lchi[nw]-this->jchi[nw]-0.5>1e-7 && this->lchi[nw]-this->jchi[nw]-0.5<1e-7)
             {
                  cout<<"Ignore ADDINFO section"<<endl;
                  this->has_so = 0;
             }
       }
       //RELBETA
       for(int nb = 0;nb < nbeta;nb++)
       {
             ifs >> this->lll[nb] >> this->jjj[nb];
             if(this->lll[nb]-this->jjj[nb]-0.5>1e-7 && this->lll[nb]-this->jjj[nb]-0.5<1e-7)
             {
                  cout<<"Ignore ADDINFO section"<<endl;
                  this->has_so = 0;
             }
       }
       return;
}


void Pseudopot_upf::print_pseudo_upf(ofstream &ofs)
{
	TITLE("Pseudopot_upf","print_pseudo_upf");
	ofs << " ==== read_pseudo_upf === " << endl;

	// print header
	ofs << " has_so: " << has_so << endl;
	ofs << " Version number : " << nv << endl;
	ofs << " Element label : " << psd << endl;
	ofs << " pp_type: " << pp_type << endl;
	ofs << " tvanp: " << tvanp << endl;
	ofs << " nlcc: " << nlcc << endl; 
	ofs << " dft: " << dft[0] << " " << dft[1] << " " << dft[2] << " " << dft[3] << endl;
	ofs << " zp: " << zp << endl;
	ofs << " etotps: " << etotps << endl;
	ofs << " ecutwfc: " << ecutwfc << endl;
	ofs << " ecutrho: " << ecutrho << endl;
	ofs << " lmax: " << lmax << endl;
	ofs << " mesh: " << mesh << endl;
	ofs << " nwfc: " << nwfc << endl;
	ofs << " nbeta: " << nbeta << endl;
	for(int i=0; i<nwfc; ++i)
	{
		ofs << " iw=" << i << " els=" << els[i] << " lchi=" << lchi[i] << " oc=" << oc[i] << endl;
	}

	ofs << " End of pseudopot_upf." << endl;

	return;

}


int Pseudopot_upf::average_p()
{
	int error = 0;
	if(!this->has_so && LSPINORB) 
	{
		error++; 
		cout<<"warning_quit! no soc upf used for lspinorb calculation, error!"<<endl; 
		return error;
	}
	//WARNING_QUIT("average_p", "no soc upf used for lspinorb calculation, error!");

	if(!this->has_so || LSPINORB) 
	{
		return error; 
	}

	int new_nbeta = 0; //calculate the new nbeta
	for(int nb=0; nb< this->nbeta; nb++)
	{
		new_nbeta++;
		if(this->lll[nb] != 0 && abs(this->jjj[nb] - this->lll[nb] - 0.5) < 1e-6) //two J = l +- 0.5 average to one
		{
			new_nbeta--;
		}
	}

	this->nbeta = new_nbeta;
	matrix dion_new;
	dion_new.create(this->nbeta, this->nbeta);

	int old_nbeta=-1;
	for(int nb=0; nb<this->nbeta; nb++)
	{
		old_nbeta++;
		int l = this->lll[old_nbeta];
		int ind=0, ind1=0;
		if(l != 0)
		{
			if(abs(this->jjj[old_nbeta] - this->lll[old_nbeta] + 0.5) < 1e-6)
			{
				if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]-0.5)>1e-6) 
					WARNING_QUIT("average_p", "error beta function 1 !");
				ind = old_nbeta +1;
				ind1 = old_nbeta;
			}
			else
			{
				if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]+0.5)>1e-6)
					WARNING_QUIT("average_p", "error beta function 2 !");
				ind = old_nbeta;
				ind1 = old_nbeta +1;
			}
			double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
			//average beta (betar)
			for(int ir = 0; ir<this->mesh;ir++)
			{
				this->beta(nb, ir) = 1.0 / (2.0 * l + 1.0) * 
						( (l + 1.0) * sqrt(this->dion(ind,ind) / vion1) *
						this->beta(ind, ir) + 
						l * sqrt(this->dion(ind1,ind1) / vion1) *
						this->beta(ind1, ir) ) ;
			}
			//average the dion matrix
			this->dion(nb, nb) = vion1;
			old_nbeta++;	
		}
		else
		{
			for(int ir = 0; ir<this->mesh;ir++)
				this->beta(nb, ir) = this->beta(old_nbeta, ir);
			this->dion(nb, nb) = this->dion(old_nbeta, old_nbeta);
		}
		this->lll[nb] = this->lll[old_nbeta]; //reset the lll index, ignore jjj index
	}

	//store the old dion and then recreate dion 
	for(int i=0;i<this->nbeta; i++)
	{
		for(int j=0;j<this->nbeta;j++)
		{
			dion_new(i,j) = this->dion(i,j);
		}
	}

	this->dion = dion_new;
//	this->dion.create(this->nbeta, this->nbeta);
//	for(int i=0;i<this->nbeta; i++)
//		for(int j=0;j<this->nbeta;j++)
//			this->dion(i,j) = dion_new(i,j);
	
	int new_nwfc = 0;
	for(int nb=0; nb<this->nwfc; nb++)
	{
		new_nwfc++;
		if(this->lchi[nb] != 0 && abs(this->jchi[nb] - this->lchi[nb] - 0.5)<1e-6)
		{
			new_nwfc--;
		}
	}

	this->nwfc = new_nwfc;
	int old_nwfc=0;
	for(int nb=0; nb<this->nwfc; nb++)
	{
		old_nwfc++;
		int l = this->lchi[old_nwfc];
		int ind=0, ind1=0;
		if(l!=0)
		{
			if(abs(this->jchi[old_nwfc] - this->lchi[old_nwfc] + 0.5) < 1e-6)
			{
				if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]-0.5)>1e-6) 
				{error++; cout<<"warning_quit! error chi function 1 !"<<endl; return error;}
//					WARNING_QUIT("average_p", "error chi function 1 !");
				ind = old_nwfc +1;
				ind1 = old_nwfc;
			}
			else
			{
				if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]+0.5)>1e-6)
				{error++; cout<<"warning_quit! error chi function 2 !"<<endl; return error;}
//					WARNING_QUIT("average_p", "error chi function 2 !");
				ind = old_nwfc;
				ind1 = old_nwfc +1;
			}
			//average chi
			for(int ir = 0; ir<this->mesh;ir++)
			{
				this->chi(nb, ir) = 1.0 / (2.0 * l + 1.0) *
					( (l+1.0)*this->chi(ind,ir) + (l*this->chi(ind1,ir)) );
				old_nwfc++;
			}
		}
		else{
			for(int ir = 0; ir<this->mesh;ir++)
				this->chi(nb, ir) = this->chi(old_nwfc, ir);
		}
		this->lchi[nb] = this->lchi[old_nwfc]; //reset lchi index
	}
	this->has_so = 0;	
	return error;
}

// Peize Lin add for bsse 2021.04.07
void Pseudopot_upf::set_empty_element(void)
{
	this->zp = 0;
	for(int ir=0; ir<mesh; ++ir)
	{
		this->vloc[ir] = 0;
	}
	for(int i=0; i<nbeta; ++i)
	{
		for(int j=0; j<nbeta; ++j)
		{
			this->dion(i,j) = 0;
		}
	}
	for(int ir=0; ir<mesh; ++ir)
	{
		this->rho_at[ir] = 0;
	}
	return;
}
