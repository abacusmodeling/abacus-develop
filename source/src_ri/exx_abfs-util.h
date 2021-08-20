#ifndef EXX_ABFS_UTIL_H
#define EXX_ABFS_UTIL_H

#include "exx_abfs.h"

#include <vector>
#include <string>

class Exx_Abfs::Util
{
public:
	static void bcast( std::vector<std::string> &v, const int rank_src, const MPI_Comm &mpi_comm );
};

#endif