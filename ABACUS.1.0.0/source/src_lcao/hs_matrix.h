#ifndef HS_MATRIX_H
#define HS_MATRIX_H

#include "../src_pw/tools.h"

// mohan add this file 2010-09-10
namespace HS_Matrix
{
	void saving_HS(const double *Hloc, const double* Sloc, bool bit, const int &out_hs);

	void save_HS(const double *H, const double *S, bool bit);
	void save_HS_complex(const complex<double> *H, const complex<double> *S, bool bit);

	void save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
		const double *nzval_H, const double *nzval_S, bool bit);
}

#endif
