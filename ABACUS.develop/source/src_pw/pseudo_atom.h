//=====================================================================
// Atomic part
// <PP_MESH>  ...  </PP_MESH>     // the mesh where the PP is defined
// <PP_NLCC>  ...  </PP_NLCC>     // non-linear core correction, optional
// <PP_RHOATOM> ...</PP_RHOATOM>  // atomic charge density
// <PP_PSWFC> ...  </PP_PSWFC>    // atomic wavfunctions
// The variables needed to describe the atoms and related quantities
// based on atom.f90
// USE parameters, ONLY : npsx, ndmx, nchix //maybe not necessary
//=====================================================================
#ifndef PSEUDO_ATOM_H
#define PSEUDO_ATOM_H

#include "tools.h"
#include "pseudopot_upf.h"
#include "pseudo_h.h"

class pseudo_atom: public pseudo_h
{
public:
	// <PP_MESH>
	double *r;        // radial logaritmic mesh, r[0:mesh-1]
	double *rab;      // derivative of the radial mesh, rab[0:mesh-1]

	//<PP_NLCC>
	double *rho_atc;  // radial core charge density, rho_atc[0:mesh-1]

	//<PP_RHOATOM>
	double *rho_at;   // radial atomic charge density, rho_at[0:mesh-1]

	// <PP_PSWFC>
	matrix chi;	  // radial atomic orbitals, chi(nchi, mesh)

	//other
	int msh;          // the point at rcut
	double rcut;      // cut-off radius

	//member functions
	pseudo_atom();
	~pseudo_atom();

	void set_pseudo_at(const Pseudopot_upf &upf);
	void print_pseudo_at(ofstream &ofs);
};

#endif
