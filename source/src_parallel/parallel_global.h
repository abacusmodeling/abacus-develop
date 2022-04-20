//==========================================================
// AUTHOR : Fang Wei, Mohan Chen
// DATE : 2008
// LAST UPDATE : 2009-3-23 mohan add GATHER_MINIMUM_DOUBLE
//==========================================================
#ifndef PARALLEL_GLOBAL_H
#define PARALLEL_GLOBAL_H

#include <complex>
using namespace std;

#ifdef __MPI
#include "mpi.h"
extern MPI_Datatype mpicomplex;
extern MPI_Op myOp;
extern MPI_Comm POOL_WORLD;
extern MPI_Comm PARAPW_WORLD;
extern MPI_Comm GRID_WORLD; //mohan add 2012-01-13
extern MPI_Comm DIAG_WORLD; //mohan add 2012-01-13
#endif

//void myProd(std::complex<double> *in,std::complex<double> *inout,int *len,MPI_Datatype *dptr);

namespace Parallel_Global
{
	//---------------------------
	// call at the very first.
	//---------------------------
	void read_mpi_parameters(int argc, char **argv);
#ifdef __MPI
	void myProd(std::complex<double> *in,std::complex<double> *inout,int *len,MPI_Datatype *dptr);
#endif
	
	//-------------------------------------------
	// call to split the "diago world"
	// the unit of first proc of each grid group 
	// us the diag world.
	// for example, if we have 64 processors,
	// and diago_proc = 4,
	// then we have 4 'grid world', which
	// have 16, 16, 16, 16 processors each,
	// and the head of each 'grid world'
	// leads to the 'diag world', diag
	// is only carried out using those 4 proc.
	//-------------------------------------------
	void split_diag_world(const int &diag_np);
	void split_grid_world(const int &diag_np);
}


#endif // GMPI
