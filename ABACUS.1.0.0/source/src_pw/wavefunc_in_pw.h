#ifndef WAVEFUNC_IN_PW_H
#define WAVEFUNC_IN_PW_H
#include "../src_pw/global.h"
#include "../src_pw/tools.h"

namespace Wavefunc_in_pw
{
	void make_table_q(std::vector<string> &orbital_files, realArray &table_local);
	void integral(const int meshr, const double *psir, const double *r,
		const double *rab, const int &l, double* table);
	
	//mohan add 2010-04-20
	double smearing(const double &energy_x,const double &ecut,const double &beta);
	void produce_local_basis_in_pw(const int &ik,ComplexMatrix &psi, const realArray &table_local); 
};
#endif
