#ifndef PDIAG_DOUBLE_H
#define PDIAG_DOUBLE_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "diag_scalapack_gvx.h"
#include "src_parallel/parallel_orbitals.h"

class Pdiag_Double : public Parallel_Orbitals
{
	public:
	Pdiag_Double();
	~Pdiag_Double();

	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	//void diago_double_begin(const int &ik, double **wfc, ModuleBase::matrix &wfc_2d,
	//	double *h_mat, double *s_mat, double *ekb);			// Peize Lin add wfc_2d 2019-01-17
	void diago_double_begin(const int &ik, ModuleBase::matrix &wfc_2d, double *h_mat, double *s_mat, double *ekb); //LiuXh add 2021-09-06, clear memory, totwfc not used now
	
	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	void diago_complex_begin(const int &ik, std::complex<double> **wfc, ModuleBase::ComplexMatrix &wfc_2d,
        std::complex<double>* ch_mat, std::complex<double>* cs_mat, double* ekb);			// Peize Lin add wfc_2d 2019-01-17

#ifdef __MPI
    void readin(const std::string& fa, const std::string& fb, const int& nlocal, double* eigen, double* eigvr);
#endif

    Diag_Scalapack_gvx diag_scalapack_gvx;			// Peize Lin add 2021.11.02


private:

#ifdef __MPI
	//void gath_eig(MPI_Comm comm,int n,double **c,double *Z);
	void gath_eig(MPI_Comm comm,int n,double *Z); //LiuXh add 2021-09-06, clear memory, totwfc not used now
	void gath_eig_complex(MPI_Comm comm,int n,std::complex<double> **c,std::complex<double> *Z, const int &ik); //mohan add 2012-01-09
	void gath_full_eig(MPI_Comm comm,int n,double **c,double *Z);
	void gath_full_eig_complex(MPI_Comm comm,int n,std::complex<double> **c, std::complex<double> *Z);
#endif

};

#endif
