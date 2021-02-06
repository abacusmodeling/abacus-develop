/*************************************************************
// Norm conserving NON-LOCAL pseudopotential (NCNL-PP)
// include pseudo_vl
// <PP_BETA>  ... </PP_BETA>       // projector
// <PP_DIJ> ... </PP_DIJ>          // d(i,j)

// Norm-conserving (atomitic) potential
// based on MODULE uspp_param
// USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx
*************************************************************/

#ifndef PSEUDO_NC_H
#define PSEUDO_NC_H

#include "tools.h"
#include "pseudopot_upf.h"
#include "pseudo_vl.h"

class pseudo_nc: public pseudo_vl
{

	public:

	// <PP_BETA>
	int *lll;       // lll(nbeta), angular momentum of the beta function
	int kkbeta;		// kkbeta(nbeta), point where the beta are zero

	// <PP_DIJ>
	matrix dion;	// dion(nbeta,nbeta)
	matrix betar;	//(nbeta, mesh), radial beta_{mu} functions

	// other
	int nh;         // number of beta functions per atomic type

	pseudo_nc();
	~pseudo_nc();

	void set_pseudo_nc(const Pseudopot_upf &upf);
	void print_pseudo_nc(ofstream &ofs);

};

#endif // PSEUDO_NC_H
