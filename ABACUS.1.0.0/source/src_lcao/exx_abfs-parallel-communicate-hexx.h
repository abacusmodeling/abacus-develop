#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_H

#include "exx_abfs-parallel.h"
#include "src_global/matrix.h"
#include "src_global/matrix_wrapper.h"
#include "src_global/complexmatrix.h"
#include "src_lcao/abfs-vector3_order.h"
#include <vector>
#include <map>
#include <deque>
#include <mpi.h>
#include <atomic>
#include <boost/dynamic_bitset.hpp>
#include <boost/mpi.hpp>

class Exx_Abfs::Parallel::Communicate::Hexx
{
public:
	void Rexx_to_Km2D(
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_exx,
		const pair<bool,bool> &io_HR_a2D );

private:
	map<size_t,map<size_t,matrix>> R_to_K(
		map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR) const;
	map<size_t,map<size_t,ComplexMatrix>> R_to_K(
		const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR,
		const size_t ik) const;
	template< typename T >
	T a2D_to_m2D( const map<size_t,map<size_t,T>> & H_a2D ) const;
	template<typename T>
	T pulay_mixing( const T &H_pulay_old, deque<T> &H_seq, const T &H_new );

	class Allreduce
	{
	public:
		Allreduce(
			const MPI_Comm &mpi_comm_in,
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local_in );
		~Allreduce();
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> exx_to_a2D();

	private:
		void set_atom_in_2D();
		void ask( const int rank_delta_now );
		void recv_data_process( const int rank_data );
		void insert_data( vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> &data_rank );
		void insert_data();
		void send_data_process( const int rank_asked );
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> get_data_local_wrapper( 
			const pair<vector<bool>,vector<bool>> &atom_asked ) const;

	private:
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> data_all;
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local;

		pair< vector<bool>, vector<bool> > atom_in_2D;

		const MPI_Comm & mpi_comm;
		int comm_sz;
		int my_rank;

		static constexpr int tag_ask  = 1;
		static constexpr int tag_data = 2;

		atomic_flag lock_insert;
		atomic<int> rank_delta;

		boost::mpi::packed_oarchive *oarp_atom_in_2D;
		vector<boost::mpi::packed_oarchive*> oarps_isend_data;
		vector<boost::mpi::packed_iarchive*> iarps_recv_data;
		vector<boost::mpi::packed_iarchive*> iarps_atom_asked;

		vector<atomic<int>*> flags_isend_data;
		vector<atomic<int>*> flags_ask_atom;
		boost::dynamic_bitset<> flags_recv_data;
	};

public:
	vector<matrix>        HK_Gamma_m2D;
	vector<ComplexMatrix> HK_K_m2D;
	vector<deque<matrix>>        HK_Gamma_m2D_pulay_seq;
	vector<deque<ComplexMatrix>> HK_K_m2D_pulay_seq;

	enum class Mixing_Mode{ No, Plain, Pulay };
	Mixing_Mode mixing_mode;
	double mixing_beta;
};

#endif