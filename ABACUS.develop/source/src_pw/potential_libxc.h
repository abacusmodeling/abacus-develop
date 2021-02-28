//==========================================================
// AUTHOR : Peize Lin
// DATE :   2017-09-14
// UPDATE : 2021-02-28
//==========================================================

#ifdef USE_LIBXC

#ifndef POTENTIAL_LIBXC_H
#define POTENTIAL_LIBXC_H

#include "src_global/matrix.h"
#include "src_global/vector3.h"
#include "src_global/global_variable.h"
#include "src_global/global_function.h"
#include <xc.h>
#include <vector>
#include <tuple>

class Potential_Libxc
{
	public:

	static void v_xc(
		const double * const * const rho_in,
		double &etxc,
		double &vtxc,
		matrix &v);
		
	private:

	static std::vector<XC(func_type)> init_func();

	// [rho, sigma, gdr] = cal_input( funcs, rho_in )
	static std::tuple< 
		std::vector<double>, 
		std::vector<double>, 
		std::vector<std::vector<Vector3<double>>> > 
	cal_input(
		const std::vector<XC(func_type)> &funcs, 
		const double * const * const rho_in );

	static int nspin0()			// 修改了soc判据后，这里可能需要修改
	{
		if     (NSPIN==1 || (NSPIN==4 && (!DOMAG && !DOMAG_Z)))		return 1;
		else if(NSPIN==2 || (NSPIN==4 && ( DOMAG ||  DOMAG_Z)))		return 2;
		else throw runtime_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}
};

#endif

#endif	// ifdef USE_LIBXC
