#ifndef POTENTIAL_LIBXC_H
#define POTENTIAL_LIBXC_H

#include "src_global/matrix.h"
#include "src_global/vector3.h"
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
	// [ rho, sigma, gdr ]
	static std::tuple< std::vector<double>, std::vector<double>, std::vector<std::vector<Vector3<double>>> > 
	cal_input( 
		const std::vector<XC(func_type)> &funcs, 
		const double * const * const rho_in );
};

#endif