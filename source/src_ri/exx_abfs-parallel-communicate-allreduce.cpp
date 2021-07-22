#include "exx_abfs-parallel-communicate-allreduce.h"

#include <thread>

Exx_Abfs::Parallel::Communicate::Allreduce::Allreduce( const MPI_Comm & mpi_comm_in )
	:mpi_comm(mpi_comm_in)
{
	MPI_Comm_size( mpi_comm, &comm_sz );
	MPI_Comm_rank( mpi_comm, &my_rank );
}


void Exx_Abfs::Parallel::Communicate::Allreduce::ask( const int rank_delta_now ) const
{
	if( rank_delta_now < comm_sz )
	{
		const int rank_ask = ( my_rank + rank_delta_now ) % comm_sz;
		MPI_Request request_ask;
		MPI_Isend( &tag_ask, 1, MPI_INT, rank_ask, tag_ask, mpi_comm, &request_ask );
	}
}