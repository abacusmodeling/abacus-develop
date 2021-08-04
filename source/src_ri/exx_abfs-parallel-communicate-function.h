#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_FUNCTION_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_FUNCTION_H

#include "exx_abfs-parallel.h"
#include <vector>
#include <utility>
#ifdef __MPI
#include <mpi.h>
#endif
using namespace std;

class Exx_Abfs::Parallel::Communicate::Function
{
public:
	static std::vector<pair<std::vector<bool>,std::vector<bool>>> get_atom_in_2D_list(const MPI_Comm &mpi_comm);
};

#endif