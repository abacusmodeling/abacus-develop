//==========================================================
// AUTHOR : Peize Lin
// DATE :   2017-09-14
// UPDATE : 2021-02-28
//==========================================================

#ifdef USE_LIBXC

#ifndef POTENTIAL_LIBXC_H
#define POTENTIAL_LIBXC_H

#include "../module_base/matrix.h"
#include "../module_base/vector3.h"
#include "../module_base/global_variable.h"
#include "../module_base/global_function.h"
#include <xc.h>
#include <vector>
#include <tuple>

class Potential_Libxc
{
	public:

	//------------------------------------------------
	// evaluate the exchange-correlation (XC) energy
	// by using the input charge density rho_in and rho_core_in
	//------------------------------------------------
	// [etxc, vtxc, v] = v_xc(...)
	static std::tuple<double,double,ModuleBase::matrix> v_xc(
		const int &nrxx, // number of real-space grid
		const int &ncxyz, // total number of charge grid
		const double &omega, // volume of cell
		const double * const * const rho_in,
		const double * const rho_core_in);
	
	static std::tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> v_xc_meta(
		const double * const * const rho_in,
		const double * const rho_core_in,
		const double * const * const kin_r_in);
		
	private:

	//-------------------------------------------
	// return the type of XC functional by 
	// calling init_func()
	//-------------------------------------------
	static std::vector<XC(func_type)> init_func();

	//------------------------------------------------
	// evaluate three quantities: rho, sigma, and gdr
	// according to the input types of XC functionals
	//------------------------------------------------
	// [rho, sigma, gdr] = cal_input(...)
	static std::tuple< 
		std::vector<double>, 
		std::vector<double>, 
		std::vector<std::vector<ModuleBase::Vector3<double>>> > 
	cal_input(
		const std::vector<XC(func_type)> &funcs, 
		const double * const * const rho_in,
		const double * const rho_core_in );

	//----------------------------
	// decide the value of spin
	//----------------------------
	static int nspin0() // may need updates from SOC
	{
		if     (GlobalV::NSPIN==1 || (GlobalV::NSPIN==4 && (!GlobalV::DOMAG && !GlobalV::DOMAG_Z)))		return 1;
		else if(GlobalV::NSPIN==2 || (GlobalV::NSPIN==4 && ( GlobalV::DOMAG ||  GlobalV::DOMAG_Z)))		return 2;
		else throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}
};

#endif

#endif	// ifdef USE_LIBXC
