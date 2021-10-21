#ifndef PDEEPKS_H
#define PDEEPKS_H

#ifdef __DEEPKS
#ifdef __MPI
#include "mpi.h"
#endif
#include "../src_pw/global.h"

//In DeePKS, we need to calculate the overlap between atomic basis and projectors {alpha}
//which is stored in matrix S_mu_alpha and calculated in build_S_descriptor
//The size of matrix if thus Nalpha * NLOCAL
//The parallelization scheme for this step is to distribute atomic orbitals

class Parallel_deepks
{
	public:
	Parallel_deepks();
	~Parallel_deepks();

	int* norb_local;
	//number of atomic basis on local processor
	int* trace_loc_orb;
	//array of size NLOCAL
	//trace_loc_orb[iw_all] = iw_local, where iw_local runs from 0 to norb_local[rank]-1
	// = -1 if orbital is not treated on current processor

	//set norb_local
	void set_nlocal(
#ifdef __MPI
		MPI_Comm MPI_WORLD
#endif
	);

	//set trace_loc_orb
	void set_loc_orb(
#ifdef __MPI
		MPI_Comm MPI_WORLD
#endif
	);

	//collect S_mu_alpha from all ranks
	void allsum_S_mu_alpha(
		int inlmax, //first dimension of S_mu_alpha
		int ndim, //secohd dimension of S_mu_alpha
		double** S_mu_alpha);

};

namespace GlobalC
{
extern Parallel_deepks ParaD;
}

#endif
#endif