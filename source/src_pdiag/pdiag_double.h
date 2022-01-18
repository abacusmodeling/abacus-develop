#ifndef PDIAG_DOUBLE_H
#define PDIAG_DOUBLE_H

#include "../src_pw/tools.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "pdiag_basic.h"
#include "diag_scalapack_gvx.h"

class Pdiag_Double : public Pdiag_Basic
{
	public:
	Pdiag_Double();
	~Pdiag_Double();

	LocalMatrix MatrixInfo;

	int nrow;
	int ncol;
	long nloc;

#ifdef __MPI
    int blacs_ctxt;
	int mpi_comm_rows, mpi_comm_cols;
    int desc[9];
#endif

	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	//void diago_double_begin(const int &ik, double **wfc, ModuleBase::matrix &wfc_2d,
	//	double *h_mat, double *s_mat, double *ekb);			// Peize Lin add wfc_2d 2019-01-17
	void diago_double_begin(const int &ik, ModuleBase::matrix &wfc_2d, double *h_mat, double *s_mat, double *ekb); //LiuXh add 2021-09-06, clear memory, totwfc not used now
	
	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	void diago_complex_begin(const int &ik, std::complex<double> **wfc, ModuleBase::ComplexMatrix &wfc_2d,
		std::complex<double> *ch_mat, std::complex<double> *cs_mat, double *ekb);			// Peize Lin add wfc_2d 2019-01-17

#ifdef __MPI
	void readin(const std::string &fa, const std::string &fb, const int &nlocal, double *eigen, double *eigvr);
#endif

	// be called
	void divide_HS_2d
	(
#ifdef __MPI
		MPI_Comm DIAG_WORLD
#endif
	);

#ifdef __MPI
	MPI_Comm comm_2D;
#endif

	int nb;
	int dim0;
	int dim1;

	Diag_Scalapack_gvx diag_scalapack_gvx;			// Peize Lin add 2021.11.02
};

#endif
