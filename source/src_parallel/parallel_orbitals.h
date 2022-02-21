#ifndef PARALLEL_ORBITALS_H
#define PARALLEL_ORBITALS_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_pdiag/pdiag_common.h"

class Parallel_Orbitals
{
    public:

    Parallel_Orbitals();
    ~Parallel_Orbitals();

    // type : Sloc(1) Hloc(2) Hloc_fixed(3)
    bool in_this_processor(const int &iw1_all, const int &iw2_all);
    
    int* trace_loc_row;
    int* trace_loc_col;
    int out_hs; // mohan add 2010-09-02
    int out_hsR; // LiuXh add 2019-07-16

    void set_trace(void);

    // use MPI_Alltoallv to convert a orbital distributed matrix to 2D-block cyclic distributed matrix
    // here are the parameters which will be used
    int sender_index_size;
    int *sender_local_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int receiver_index_size;
    int *receiver_global_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;
    
    LocalMatrix MatrixInfo;

	int nrow;
	int ncol;
	long nloc;

#ifdef __MPI
    int blacs_ctxt;
	int mpi_comm_rows, mpi_comm_cols;
    int desc[9];
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

    // mohan add 2010-09-10
	// output local wave functions.
	// put it here because if we 
	// use HPSEPS, the wave functions
	// is needed to be collected first.
    int out_lowf;
    double** Z_LOC; //xiaohui add 2014-06-19
    
protected: //also used in diago
    int testpb;
    bool alloc_Z_LOC; //xiaohui add 2014-12-22
    int lastband_in_proc;
	int lastband_number; 
	int* loc_sizes;
    int loc_size;

    void set_parameters(void);

#ifdef __MPI
	void mpi_creat_cart(MPI_Comm *comm_2D, int prow, int pcol);

	void mat_2d(MPI_Comm vu, const int &M_A, const int &N_A, const int & NB, LocalMatrix &loc_A);

	void data_distribution(
		MPI_Comm comm_2D, 
		const std::string &file,
		const int &n,
		const int &NB,
		double *A,
		const LocalMatrix &loc_A);

#endif

};

#endif
