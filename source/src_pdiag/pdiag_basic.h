#ifndef PDIAG_BASIC_H
#define PDIAG_BASIC_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

class Pdiag_Basic
{
	public:
	Pdiag_Basic();
	~Pdiag_Basic();

	// mohan add 2010-09-10
	// output local wave functions.
	// put it here because if we 
	// use HPSEPS, the wave functions
	// is needed to be collected first.
	int out_lowf;
	double** Z_LOC; //xiaohui add 2014-06-19

	protected:

	void set_parameters(void);

#ifdef __MPI
	void mpi_creat_cart(MPI_Comm *comm_2D, int prow, int pcol);

	void mat_2d(MPI_Comm vu, const int &M_A, const int &N_A, const int & NB, LocalMatrix &loc_A);

	void data_distribution(
		MPI_Comm comm_2D, 
		const string &file,
		const int &n,
		const int &NB,
		double *A,
		const LocalMatrix &loc_A);

	void gath_eig(MPI_Comm comm,int n,double **c,double *Z);
	void gath_eig_complex(MPI_Comm comm,int n,complex<double> **c,complex<double> *Z, const int &ik); //mohan add 2012-01-09
	void gath_full_eig(MPI_Comm comm,int n,double **c,double *Z);
	void gath_full_eig_complex(MPI_Comm comm,int n,complex<double> **c, complex<double> *Z);
#endif

	int lastband_in_proc;
	int lastband_number; 
	int* loc_sizes;
	int loc_size;
	
	protected:

	int testpb;
	bool alloc_Z_LOC; //xiaohui add 2014-12-22


};

#endif
