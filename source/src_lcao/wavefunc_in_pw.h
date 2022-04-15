#ifndef WAVEFUNC_IN_PW_H
#define WAVEFUNC_IN_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"
#include "../module_base/realarray.h"
#include "../module_base/vector3.h"

//---------------------------------------------------
// FUNCTION: expand the local basis sets into plane
// wave basis sets
//---------------------------------------------------
namespace Wavefunc_in_pw
{

	void make_table_q(
		std::vector<std::string> &orbital_files, 
		ModuleBase::realArray &table_local);

	void integral(
		const int meshr, // number of mesh points 
		const double *psir,
		const double *r,
		const double *rab, 
		const int &l, 
		double* table);
	
	//mohan add 2010-04-20
	double smearing(
		const double &energy_x,
		const double &ecut,
		const double &beta);

	void produce_local_basis_in_pw(
		const int &ik,
		ModuleBase::ComplexMatrix &psi, 
		const ModuleBase::realArray &table_local);

 	void produce_local_basis_q_in_pw(
		const int &ik,
		ModuleBase::ComplexMatrix &psi, 
		const ModuleBase::realArray &table_local, 
		ModuleBase::Vector3<double> q);     // pengfei 2016-11-23

}
#endif
