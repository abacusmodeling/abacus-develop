#ifndef EXX_ABFS_IO_TEMPLATE_H
#define EXX_ABFS_IO_TEMPLATE_H

#include "src_lcao/exx_abfs-io.h"
#include "src_global/boost_serialization.h"
#include <string>
#include <fstream>

template<typename T>
void Exx_Abfs::IO::output_binary( const T &data, const std::string &file_name )
{
	std::ofstream ofs(file_name,std::ios::binary);
	boost::archive::binary_oarchive oa(ofs);
	oa << data;
	ofs.close();
}

template<typename T>
T Exx_Abfs::IO::input_binary( const std::string &file_name )
{
	std::ifstream ifs(file_name,std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	T data;
	ia >> data;
	ifs.close();
	return data;
}

template<typename T>
void Exx_Abfs::IO::output_text( const T &data, const std::string &file_name )
{
	std::ofstream ofs(file_name);
	boost::archive::text_oarchive oa(ofs);
	oa << std::move(data);
	ofs.close();
}

template<typename T>
T Exx_Abfs::IO::input_text( const std::string &file_name )
{
	std::ifstream ifs(file_name);
	boost::archive::text_iarchive ia(ifs);
	T data;
	ia >> data;
	ifs.close();
	return data;
}

template<typename T>
void Exx_Abfs::IO::bcast( T &data, const int rank_src, MPI_Comm mpi_comm )
{
	int my_rank;	MPI_Comm_rank( mpi_comm, &my_rank );
	if(MY_RANK==rank_src)
	{
		boost::mpi::packed_oarchive oar(mpi_comm);
		oar << data;
		const int data_size = oar.size();
		MPI_Bcast( const_cast<int*>(&data_size), 1, MPI_INT, rank_src, mpi_comm );
		MPI_Bcast( const_cast<void*>(oar.address()), oar.size(), MPI_PACKED, rank_src, mpi_comm );
	}
	else
	{
		int data_size;
		MPI_Bcast( &data_size, 1, MPI_INT, rank_src, mpi_comm );
		boost::mpi::packed_iarchive iar(mpi_comm);
		iar.resize(data_size);
		MPI_Bcast( iar.address(), data_size, MPI_PACKED, rank_src, mpi_comm );
		iar >> data;
	}
}

#endif