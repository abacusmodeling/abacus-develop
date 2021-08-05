#ifndef EXX_ABFS_PARALLEL_DISTRIBUTE_ORDER_H
#define EXX_ABFS_PARALLEL_DISTRIBUTE_ORDER_H

#include "exx_abfs-parallel.h"

#include "abfs.h"

#include <vector>
#include <utility>

class Exx_Abfs::Parallel::Distribute::Order
{
public:
	static std::vector<std::pair<size_t,size_t>> distribute(
		const double rmesh_times );
};

#endif