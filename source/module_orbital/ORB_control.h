#ifndef ORB_CONTROL_H 
#define ORB_CONTROL_H 

#include "src_parallel/parallel_orbitals.h"
#include "src_parallel/parallel_global.h"
#include "ORB_gen_tables.h"
#include "ORB_read.h"

class ORB_control 
{
	public:

    ORB_control();
    ~ORB_control();

	//first step: read orbital file
	void read_orb_first(
	std::ofstream &ofs_in, 
	LCAO_Orbitals &orb,
	const int &ntype, // mohan add 2021-04-26
	const int &lmax, // mohan add 2021-04-26 
	const double &lcao_ecut_in, // mohan add 2021-04-16
	const double &lcao_dk_in, // mohan add 2021-04-16
	const double &lcao_dr_in, // mohan add 2021-04-16
	const double &lcao_rmax_in, // mohan add 2021-04-16
	const int &out_descriptor,
	const int &out_r_matrix,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank // mohan add 2021-04-26
	);
    /// Generate the S(overlap),T,NL matrix.
    void set_orb_tables(
		std::ofstream &ofs_in,
		ORB_gen_tables &OGT, 
		LCAO_Orbitals &orb,
		const double &lat0,
		const int &out_descriptor,
		const int &Lmax_exx,
		const int &nprojmax, 
		const int* nproj,
		const Numerical_Nonlocal* beta_); 

    void clear_after_ions(
		ORB_gen_tables &OGT, 
		LCAO_Orbitals &orb,
		const int &out_descriptor,
        const int* nproj_);

    void setup_2d_division(void);
    
#ifdef __MPI
    void readin(const std::string& fa, const std::string& fb, const int& nlocal, double* eigen, double* eigvr);
#endif

    Parallel_Orbitals ParaV;

private:

    void divide_HS_2d
	(
#ifdef __MPI
		MPI_Comm DIAG_WORLD
#endif
    );
    
#ifdef __MPI
	int mpi_comm_rows, mpi_comm_cols;
#endif
    
    void set_parameters(void);
    
    void set_trace(void);

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
