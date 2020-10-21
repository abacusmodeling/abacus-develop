/*****************************************************************
// General information part
// base on pseudo_type.f90,
// <PP_INFO>  ...  </PP_INFO>     // for reader's information only
// <PP_HEADER> ... </PP_HEADER>   // general information about the pseudopotential
*******************************************************************/
#ifndef PSEUDOH_H
#define PSEUDOH_H

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include "pseudopot_upf.h"

class pseudo_h
{
public:

	//<PP_HEADER>
	bool has_so;		// if .true. includes spin-orbit
	int  nv;			// UPF file version number
	string psd;			// Element label
	string pp_type;		// Pseudo type ( NC or US )
	bool tvanp;			// .true. if Ultrasoft
	bool nlcc;			// Non linear core corrections(bool)
	string dft[4];		// Exch-Corr type
	int  zv;			// z valence
	double etotps;		// total energy
	double ecutwfc;		// suggested cut-off for wfc
	double ecutrho;		// suggested cut-off for rho
	int lmax;			// maximum angular momentum component
	int mesh;			// number of point in the radial mesh
	int nchi;			// nwfc,number of wavefunctions
	int nbeta;			// number of projectors
	string *els;		// els[nchi]
	int *lchi;			// lchi[nchi]
	double *oc;			// oc[nchi]

	/*  <pp_info> */
//  char generated[80];	  // author
//  char date_author[80]; // Misc info
//  char comment[80];	  //
//  int * nn;		  // nn(nwfc)
//  double * rcut;	  // cut-off radius(nwfc)	??
//  double * rcutus;	  // cut-off ultrasoft radius (nwfc)
//  double * epseu;	  // energy (nwfc)
//  double  xmin;         // the minimum x of the linear mesh
//  double  rmax;         // the maximum radius of the mesh
//  double  zmesh;        // the nuclear charge used for mesh
//  double  dx;           // the deltax of the linear mesh
//  double *jchi;			// jchi(nwfc)
	double *jjj;	  // total angual momentum, jjj[nbeta]
	double *jchi;	//jchi(nwfc), added by zhengdy-soc
	int *nn;

	pseudo_h();
	~pseudo_h();

	// member functions
	void set_pseudo_h(const Pseudopot_upf &upf);
	void print_pseudo_h(ofstream &ofs);
};

#endif // PSEUDOH_H

/*
// pseudo-potential that has spin-orboit term
// We are not going to implement this PP in the first version
// USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx

#ifndef PSEUDO_SO_H
#define PSEUDO_SO_H

class pseudo_so:public pseudo_us
{
	double *jjj;	//(nbrx,npsx), &! total angular momentum of the beta function
//	pseudo_us uspp;		// if SO+USPP

	pseudo_so();
	~pseudo_so();

	void set_pseudo_upf();

};
#endif // PSEUDO_SO_H
*/
