/************************************************************
// Ultrasoft pseudopotential (US-PP)
// include pseudo_nc
// <PP_QIJ>  ...  </PP_QIJ>     //psudo-charge

//We are not going to implement this PP in the first version
//but we still construct the PP here for future reference

// based on MODULE uspp_param
// USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx
#ifndef PSEUDO_US_H
#define PSEUDO_US_H
************************************************************/

#ifndef PSEUDO_US_H
#define PSEUDO_US_H

#include "tools.h"
#include "pseudopot_upf.h"
#include "pseudo_nc.h"

class pseudo_us: public pseudo_nc
{
public:
	// <PP_QIJ>
	int nqf;
	//  <PP_RINNER>
	double *rinner;		// rinner(0:2*lmax)
	matrix qqq;			// qqq(nbeta,nbeta)
	realArray qfunc;	// qfunc(mesh,nbeta,nbeta)

	// <PP_QFCOEF>
	realArray qfcoef;	// qfcoef(nqf,0:2*lmax,nbeta,nbeta)

	// other ?
	int nqlc;			// number of angular momenta in Q

	pseudo_us();
	~pseudo_us();

	void set_pseudo_us(const Pseudopot_upf &upf);
	void print_pseudo_us(ofstream &ofs);

};

#endif //PSEUDO_US_H
