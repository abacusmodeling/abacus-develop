//==========================================================
// AUTHOR : mohan, Lixin He
// DATE : 2009-12-16
//==========================================================
#ifndef GGA_PW_H 
#define GGA_PW_H

#include "tools.h"
#include "../src_parallel/parallel_global.h"
//==========================================================
// Calculate the charge gradient using plane wave basis. 
//==========================================================
namespace GGA_PW
{
	void gradcorr(double &etxc, double &vtxc, matrix &v);
	void grad_rho( const complex<double> *rhog, double **gdr );
	void grad_dot(double **h, double *dh);

}

#endif
