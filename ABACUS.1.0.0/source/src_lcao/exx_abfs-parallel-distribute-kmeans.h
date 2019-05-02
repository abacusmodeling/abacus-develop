#ifndef EXX_ABFS_PARALLEL_DISTRIBUTE_KMEANS_H
#define EXX_ABFS_PARALLEL_DISTRIBUTE_KMEANS_H

#include "exx_abfs-parallel.h"
#include "src_global/vector3.h"
#include <vector>
#include <utility>
#include <mpi.h>

class Exx_Abfs::Parallel::Distribute::Kmeans
{
public:

	static vector<pair<size_t,size_t>> distribute( const MPI_Comm & mpi_comm, const int multiple_core=1 );

private:

	struct Cluster
	{
		Vector3<double> tau;
		size_t size;
		Vector3<double> tau_sum;
	};

	struct Atom
	{
		Vector3<double> tau;
		int center;
		double distance;
	};
	
	static vector<Atom> cluster( const int Nc );
};

#endif