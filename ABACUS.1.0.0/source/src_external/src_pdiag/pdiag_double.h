#ifndef PDIAG_DOUBLE_H
#define PDIAG_DOUBLE_H

#include "../../src_pw/tools.h"
#include "pdiag_basic.h"

class Pdiag_Double : public Pdiag_Basic
{
	public:
	Pdiag_Double();
	~Pdiag_Double();

	LocalMatrix MatrixInfo;

	int nrow;
	int ncol;
	int nloc;

	void diago_double_begin(const int &ik, double **c, 
		double *h_mat, double *s_mat, double *ekb);
	
	void diago_complex_begin(const int &ik, complex<double> **c,
		complex<double> *ch_mat, complex<double> *cs_mat, double *ekb);

#ifdef __MPI
	void readin(const string &fa, const string &fb, const int &nlocal, double *eigen, double *eigvr);
#endif

	void divide_HS_2d
	(
#ifdef __MPI
		MPI_Comm DIAG_WORLD
#endif
	);

	private:

#ifdef __MPI
	MPI_Comm comm_2D;
#endif

	int nb;
	int dim0;
	int dim1;

};

#endif
