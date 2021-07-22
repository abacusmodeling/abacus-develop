#ifndef PSEUDOPOT_UPF_H
#define PSEUDOPOT_UPF_H

#include <string>
#include "../module_base/matrix.h"
#include "../src_pw/tools.h"
#include "../src_io/output.h"

using namespace std;

class Pseudopot_upf
{
public:
	//PP_INFO
	//PP_HEADER
	//PP_MESH
	//PP_NLCC
	//PP_LOCAL
	//PP_NONLOCAL
	//PP_PSWFC
	//PP_PSRHOATOM
	//addinfo

	Pseudopot_upf();
	~Pseudopot_upf();

	bool has_so;        // if .true. includes spin-orbit

	int nv;             // header_1 // UPF file version number
	string psd;			// header_2 // Element label
	string pp_type;		// header_3 // Pseudo type ( NC or US )
	bool tvanp;         // header_4 // .true. if Ultrasoft 
	bool nlcc;		    // header_5 // Non linear core corrections
	string dft[4];		// header_6 // Exch-Corr type
	int  zp;            // header_7 // z valence 
	double etotps;		// header_8 // total energy
	double ecutwfc;		// header_9 // suggested cut-off for wfc
	double ecutrho;		// header_10 // suggested cut-off for rho
	int lmax;			// header_11 // maximum angular momentum component
	int mesh;			// header_12 // number of point in the radial mesh
	int nwfc;			// header_13 // number of wavefunctions
	int nbeta;			// header_14 // number of projectors
	string *els;		// header_15 // els(nwfc):label for the i-th atomic orbital (4s, 4p, etc)
	int *lchi;			// header_16 // lchi(nwfc):angular momentum
	double *oc;			// header_17 // oc(nwfc)

	// need 'new' and 'delete'
	double *r;			// mesh_1 // r(mesh)
	double *rab;		// mesh_2 // rab(mesh)
	double *rho_atc;	// nlcc_1 // rho_atc(mesh) "Nonlinear Core Correction"
	double *vloc;		// local_1 // vloc(mesh)
	matrix chi;			// pswfc_1 // chi(nwfc,mesh)
	double *rho_at;		// psrhoatom_1 // rho_at(mesh)

	int *lll;       	// nl_1  // lll(nbeta):angular momentum of projector i
	int *kkbeta;    	// nl_2  // kkbeta(nbeta):number of mesh points for projector i (must be .le.mesh )
	matrix beta;		// nl_3  // beta(nbeta,mesh)
	matrix dion;		// nl_4  //dion(nbeta,nbeta)
						         //the D_ij factors (Ry^{-1}) of the nonlocal PP:
						         //V_NL = \sum_{i,j} D_{i,j} |\beta_i><\beta_j|
	//added by zhengdy-soc
	int *nn;                // nn(nwfc) quantum number of wfc
	double *jchi;           // jchi(nwfc) j=l+1/2 or l-1/2 of wfc
	double *jjj;            // jjj(nbeta) j=l+1/2 or l-1/2 of beta

	int  nd; 			// nl_5 // Number of nonzero Dij
	
	// the followings are for the vwr format
	int spd_loc; 
	int iTB_s;
	int iTB_p;
	int iTB_d;
	double* vs;// local pseudopotential for s, unit is Hartree 
	double* vp;// local pseudopotential for p 
	double* vd;// local pseudopotential for d
	double* ws;// wave function for s
	double* wp;// wave function for p
	double* wd;// wave function for d

	// return error
	int init_pseudo_reader(const string &fn);
	void print_pseudo_upf(ofstream &ofs);

	bool functional_error;//xiaohui add 2015-03-24
	int average_p(); //zhengdy add 2020-10-20
	void set_empty_element();		// Peize Lin add for bsse 2022.04.07

private:

	int set_pseudo_type(const string &fn);
	string& trim(string &in_str);
	string  trimend(string &in_str);

	int read_pseudo_upf(ifstream &ifs);
	int read_pseudo_vwr(ifstream &ifs);
	int read_pseudo_upf201(ifstream &ifs);
	void read_pseudo_header(ifstream &ifs);
	void read_pseudo_mesh(ifstream &ifs);
	void read_pseudo_nlcc(ifstream &ifs);
	void read_pseudo_local(ifstream &ifs);
	void read_pseudo_nl(ifstream &ifs);
	void read_pseudo_pswfc(ifstream &ifs);
	void read_pseudo_rhoatom(ifstream &ifs);
	void read_pseudo_addinfo(ifstream &ifs);
	void read_pseudo_so(ifstream &ifs);
	//string get_string( char ss[]);
	//int get_int( char ss[]);
	//double get_double( char ss[]);
	//void get_char( string ss);
	void getnameval(ifstream&, int&, string * , string *);
};

#endif //pseudopot_upf class
