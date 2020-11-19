/**********************************************************************
// LOCAL pseudopotential (PP)
// include pseudo_info,pseudo_atom
// <PP_LOCAL> ...  </PP_LOCAL>    // local part of the potential

// Local part of the (atomitic) potential
// based on MODULE uspp_param
// USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx

**********************************************************************/

#ifndef PSEUDO_VL_H
#define PSEUDO_VL_H

#include "pseudopot_upf.h"
#include "pseudo_atom.h"

class pseudo_vl: public pseudo_atom
{
public:
	double *vloc_at;	// [mesh], local potential( = pseudopot_upf.vloc )

	pseudo_vl();
	~pseudo_vl();

	void set_pseudo_vl(const Pseudopot_upf &upf);
	void print_pseudo_vl(ofstream &ofs);
};

#endif // PSEUDO_VL_H
