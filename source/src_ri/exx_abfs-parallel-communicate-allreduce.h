#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_ALLREDUCE_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_ALLREDUCE_H

#include "exx_abfs-parallel.h"
#include <functional>
#include <boost/mpi.hpp>
#include <atomic>

/*
	Here are 3 kinds of parallel:
	1. exx
	2. 2D for diag
	3. grid
*/
class Exx_Abfs::Parallel::Communicate::Allreduce
{
public:

	Allreduce( const MPI_Comm & mpi_comm_in );
	
	template< typename T >
	T allreduce(
		T & data_local,
		std::function<void(T&,T&)> & insert_function);			// insert_function( data_all, data_rank )


private:

	void ask( const int rank_delta_now ) const;

	template< typename T >
	void allreduce_thread(
		boost::mpi::packed_iarchive &iar, 
		std::function<void(T&,T&)> & insert_function, 
		T &data_all,
		atomic<int> &rank_delta, 
		atomic_flag &insert_lock );
		
	int comm_sz;
	int my_rank;
	const MPI_Comm & mpi_comm;
	
	static constexpr int tag_ask  = 1;
	static constexpr int tag_data = 2;
};


#endif