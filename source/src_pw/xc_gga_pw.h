#ifndef GGA_PW_H 
#define GGA_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/vector3.h"
#include "../src_parallel/parallel_global.h"

//==========================================================
// Calculate the charge gradient using plane wave basis. 
// UPDATE : Peize Lin change all double** to ModuleBase::Vector3<double>* 2017-12-20
//==========================================================

namespace GGA_PW
{
	void gradcorr(double &etxc, double &vtxc, ModuleBase::matrix &v);
	void grad_rho( const std::complex<double> *rhog, ModuleBase::Vector3<double> *gdr );
	void grad_wfc( const std::complex<double> *rhog, const int ik, std::complex<double> **grad, const int npw );
	void grad_dot( const ModuleBase::Vector3<double> *h, double *dh);
	void noncolin_rho(double *rhoout1,double *rhoout2,double *seg);
}

#endif
