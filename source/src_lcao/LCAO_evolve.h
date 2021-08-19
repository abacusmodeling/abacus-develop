#ifndef EVOLVE_LCAO_MATRIX_H
#define EVOLVE_LCAO_MATRIX_H

#include "../src_pw/tools.h"

class Evolve_LCAO_Matrix
{
	public:
	Evolve_LCAO_Matrix();
	~Evolve_LCAO_Matrix();
	
	void evolve_complex_matrix(const int &ik, std::complex<double> **c, std::complex<double> **c_init)const;

	private:

	void using_LAPACK_complex(const int &ik, std::complex<double> **c, std::complex<double> **c_init)const;
#ifdef __MPI
	int using_ScaLAPACK_complex(const int &ik, std::complex<double>** c, std::complex<double>** c_init)const;
#endif
};

#endif
