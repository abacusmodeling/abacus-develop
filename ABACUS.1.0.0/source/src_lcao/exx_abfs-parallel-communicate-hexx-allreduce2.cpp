#include "exx_abfs-parallel-communicate-hexx.h"
#include "src_global/boost_serialization.h"

#include "src_pw/global.h"
#include "src_lcao/global_fp.h"

#include <mpi.h>
#include <thread>


#include "src_external/src_test/test_function.h"

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::init(
	const MPI_Comm &mpi_comm_in,
	const set<pair<size_t,size_t>> &H_atom_pairs_core)
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::init");
	
	mpi_comm = mpi_comm_in;
	MPI_Comm_size(mpi_comm, &comm_sz);
	MPI_Comm_rank(mpi_comm, &my_rank);
	
	atom_in_2D_list = get_atom_in_2D_list();
	send_size_list = get_send_size_list(H_atom_pairs_core, atom_in_2D_list);
	recv_size = ParaO.nrow * ParaO.ncol * sizeof(double);
}

vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> 
Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::exx_to_a2D(
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local)
{
ofstream ofs("Hmpi_"+TO_STRING(my_rank));

	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::exx_to_a2D");
	
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> data_all(NSPIN);
ofs<<__LINE__<<endl;

	vector<atomic<Flag_Send>> flags_send(comm_sz);
	vector<atomic<Flag_Recv>> flags_recv(comm_sz);
	init_elec(flags_send, flags_recv);
ofs<<__LINE__<<endl;

	vector<boost::mpi::packed_oarchive*> oarps_isend(comm_sz,nullptr);
	vector<boost::mpi::packed_iarchive*> iarps_irecv(comm_sz,nullptr);
ofs<<__LINE__<<endl;

	int rank_send_now = my_rank;
	auto rank_send_next = [&]()->int{ return (rank_send_now+1)%comm_sz; };
ofs<<__LINE__<<endl;

	atomic_flag lock_insert = ATOMIC_FLAG_INIT;
	vector<thread> threads;
	vector<MPI_Request>requests_isend(comm_sz);
	vector<MPI_Request>requests_irecv(comm_sz);
ofs<<__LINE__<<endl;

	while(!finish_judge(flags_send, flags_recv))
	{
ofs<<__LINE__<<endl;
		if(rank_send_next()!=my_rank && memory_enough(rank_send_next(), flags_send))
		{
			rank_send_now = rank_send_next();
ofs<<__LINE__<<"\t"<<rank_send_now<<endl;
			flags_send[rank_send_now] = Flag_Send::begin_oar;
			threads.push_back(std::thread(
				&Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::send_data_process, this,
				rank_send_now, std::cref(data_local), std::ref(oarps_isend), std::ref(flags_send) ));
		}
		for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		{
			if( flags_send[rank_send] == Flag_Send::finish_oar )
			{
ofs<<__LINE__<<"\t"<<rank_send<<endl;
				if(MPI_Isend( oarps_isend[rank_send]->address(), oarps_isend[rank_send]->size(), MPI_PACKED, rank_send, 0, mpi_comm, &requests_isend[rank_send] )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				flags_send[rank_send] = Flag_Send::begin_isend;
			}
		}
		for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		{
			if( flags_send[rank_send] == Flag_Send::begin_isend )
			{
ofs<<__LINE__<<"\t"<<rank_send<<endl;
				int flag_finish;
				if(MPI_Test( &requests_isend[rank_send], &flag_finish, MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				if(flag_finish)
				{
					delete oarps_isend[rank_send];	oarps_isend[rank_send] = nullptr;
					flags_send[rank_send] = Flag_Send::finish_isend;
				}
			}
		}
ofs<<__LINE__<<endl;

		{
			MPI_Status status;
			int flag_message;
			if(MPI_Iprobe( MPI_ANY_SOURCE, 0, mpi_comm, &flag_message, &status )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			if(flag_message)
			{
				int message_size;
				if(MPI_Get_count( &status, MPI_PACKED, &message_size )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				const int rank_recv = status.MPI_SOURCE;
ofs<<__LINE__<<"\t"<<rank_recv<<endl;
				iarps_irecv[rank_recv] = new boost::mpi::packed_iarchive{mpi_comm};
				iarps_irecv[rank_recv]->resize(message_size);
				if(MPI_Irecv( iarps_irecv[rank_recv]->address(), message_size, MPI_PACKED, rank_recv, 0, mpi_comm, &requests_irecv[rank_recv] )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				flags_recv[rank_recv] = Flag_Recv::begin_irecv;
			}
		}
		for(int rank_recv=0; rank_recv!=comm_sz; ++rank_recv)
		{
			if(flags_recv[rank_recv] == Flag_Recv::begin_irecv)
			{
				int flag_finish;
				if(MPI_Test( &requests_irecv[rank_recv], &flag_finish, MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				if(flag_finish)
				{
ofs<<__LINE__<<"\t"<<rank_recv<<endl;
					flags_recv[rank_recv] = Flag_Recv::begin_iar;
					threads.push_back(std::thread(
						&Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::recv_data_process, this,
						rank_recv, std::ref(data_all), std::ref(iarps_irecv), std::ref(flags_recv), std::ref(lock_insert) ));
				}
			}
		}
	}
ofs<<__LINE__<<endl;

	while( lock_insert.test_and_set() );
	insert_data(data_local, data_all);
	lock_insert.clear();
ofs<<__LINE__<<endl;
	
	for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
	{
		if( flags_send[rank_send] == Flag_Send::begin_isend )
		{
			if(MPI_Wait( &requests_isend[rank_send], MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			delete oarps_isend[rank_send];	oarps_isend[rank_send] = nullptr;
			flags_send[rank_send] = Flag_Send::finish_isend;
		}
	}
	
	for(std::thread &t : threads)
		t.join();	
	
	return data_all;
}

vector<pair<vector<bool>,vector<bool>>>
Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_atom_in_2D_list() const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_atom_in_2D_list");
	
	bool atom_in_2D[ucell.nat*2];
	for(int i=0; i<ucell.nat*2; ++i)
		atom_in_2D[i]=false;
	for(int iwt=0; iwt<NLOCAL; ++iwt)
	{
		const int iat = ucell.iwt2iat[iwt];
		if( ParaO.trace_loc_row[iwt]>=0 )
			atom_in_2D[iat] = true;
		if( ParaO.trace_loc_col[iwt]>=0 )
			atom_in_2D[ucell.nat+iat] = true;
	}
	
	bool atom_in_2D_list_tmp[ucell.nat*2*comm_sz];
	if(MPI_Allgather( atom_in_2D, ucell.nat*2, MPI_BYTE, atom_in_2D_list_tmp, ucell.nat*2, MPI_BYTE, mpi_comm )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));

	vector<pair<vector<bool>,vector<bool>>> atom_in_2D_list(comm_sz, {vector<bool>(ucell.nat),vector<bool>(ucell.nat)});
	for(int rank=0; rank<comm_sz; ++rank)
		for(int iat=0; iat<ucell.nat; ++iat)
		{
			atom_in_2D_list[rank].first[iat] = atom_in_2D_list_tmp[rank*ucell.nat*2+iat];
			atom_in_2D_list[rank].second[iat] = atom_in_2D_list_tmp[rank*ucell.nat*2+ucell.nat+iat];
		}	
	return atom_in_2D_list;
}

// the upper limit of size sended to each process
vector<size_t> Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_send_size_list(
	const set<pair<size_t,size_t>> &H_atom_pairs_core,
	const vector<pair<vector<bool>,vector<bool>>> &atom_in_2D_list) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_send_size_list");
	vector<size_t> send_size_list(comm_sz,0);
	for(const auto &H_atom_pair : H_atom_pairs_core)
	{
		const size_t iat1 = H_atom_pair.first;
		const size_t iat2 = H_atom_pair.second;
		for(int rank_send=0; rank_send<comm_sz; ++rank_send)
		{
			if(atom_in_2D_list[rank_send].first[iat1] && atom_in_2D_list[rank_send].second[iat2])
				send_size_list[rank_send] += ucell.atoms[ucell.iat2it[iat1]].nw * ucell.atoms[ucell.iat2it[iat2]].nw * sizeof(double);
		}
	}
	return send_size_list;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::init_elec(
	vector<atomic<Flag_Send>> &flags_send,
	vector<atomic<Flag_Recv>> &flags_recv)
{
	for(atomic<Flag_Send> &flag_send : flags_send)
		flag_send = Flag_Send::undo;
	flags_send[my_rank] = Flag_Send::finish_isend;
	
	for(atomic<Flag_Recv> &flag_recv : flags_recv)
		flag_recv = Flag_Recv::undo;
	flags_recv[my_rank] = Flag_Recv::finish_iar;
}

bool Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::finish_judge(
	const vector<atomic<Flag_Send>> &flags_send,
	const vector<atomic<Flag_Recv>> &flags_recv) const
{
	for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		if((flags_send[rank_send]!=Flag_Send::begin_isend) && (flags_send[rank_send]!=Flag_Send::finish_isend))
			return false;
	for(int rank_recv=0; rank_recv!=comm_sz; ++rank_recv)
		if(flags_recv[rank_recv]!=Flag_Recv::finish_iar)
			return false;
	return true;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::send_data_process(
	const int rank_send_now,
	const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local,
	vector<boost::mpi::packed_oarchive*> &oarps_isend,
	vector<atomic<Flag_Send>> &flags_send)
{
	const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>>
		matrix_wrapped = get_data_local_wrapper(rank_send_now, data_local);
	oarps_isend[rank_send_now] = new boost::mpi::packed_oarchive(MPI_COMM_WORLD);
	*oarps_isend[rank_send_now] << matrix_wrapped;
	flags_send[rank_send_now] = Flag_Send::finish_oar;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::recv_data_process(
	const int rank_recv,
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all,
	vector<boost::mpi::packed_iarchive*> &iarps_irecv,
	vector<atomic<Flag_Recv>> &flags_recv,
	atomic_flag &lock_insert)
{
	auto vector_empty = []( const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> & v ) -> bool
	{
		for( const auto &i : v )
			if(!i.empty())	return false;
		return true;
	};

	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> data_rank;
	*iarps_irecv[rank_recv] >> data_rank;
	delete iarps_irecv[rank_recv];	iarps_irecv[rank_recv]=nullptr;
	flags_recv[rank_recv] = Flag_Recv::finish_iar;

	if(!vector_empty(data_rank))
	{
		while( lock_insert.test_and_set() );
		insert_data(data_rank, data_all);
		lock_insert.clear();
	}

	#ifdef MATRIX_WRAPPER_TIANHE2
	for( auto &data_rank_is : data_rank )
		for( auto &data_rank_A : data_rank_is )
			for( auto &data_rank_B : data_rank_A.second )
				for( auto &data_rank_C : data_rank_B.second )
					if(!data_rank_C.second.c)
						delete[] data_rank_C.second.c;
	#endif
}


bool Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::memory_enough(
	const int rank_send_next,
	const vector<atomic<Flag_Send>> &flags_send) const
{
	size_t memory_need = recv_size;
	for(int rank_send=0; rank_send<comm_sz; ++rank_send)
		if(flags_send[rank_send]==Flag_Send::begin_oar)
			memory_need += send_size_list[rank_send];
	const size_t memory_available = MemAvailable()*1024;
	return (memory_available-memory_need>send_size_list[rank_send_next]) ? true : false;
}

vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>>
Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_data_local_wrapper(
	const int rank_send_now,
	const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local) const
{
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> mw(NSPIN);
	for( int is=0; is!=NSPIN; ++is )
	{
		auto &data_local_is = data_local[is];
		auto &mw_is = mw[is];
		for( auto &data_local_A : data_local_is )
		{
			const size_t iat1 = data_local_A.first;
			if( !atom_in_2D_list[rank_send_now].first[iat1] )	continue;
			auto &mw_A = mw_is[iat1];
			for( auto &data_local_B : data_local_A.second )
			{
				const size_t iat2 = data_local_B.first;
				if( !atom_in_2D_list[rank_send_now].second[iat2] )	continue;
				auto &mw_B = mw_A[iat2];
				for( auto &data_local_C : data_local_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_local_C.first;
					mw_B[box2] = data_local_C.second;
				}
			}
		}
	}
	return mw;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::insert_data(
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local,
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all)
{
	for( int is=0; is!=NSPIN; ++is )
	{
		auto &data_local_is = data_local[is];
		auto &data_all_is = data_all[is];
		for( auto &data_local_A : data_local_is )
		{
			const size_t iat1 = data_local_A.first;
			if( !atom_in_2D_list[my_rank].first[iat1] )	continue;
			auto &data_all_A = data_all_is[iat1];
			for( auto &data_local_B : data_local_A.second )
			{
				const size_t iat2 = data_local_B.first;
				if( !atom_in_2D_list[my_rank].second[iat2] )	continue;
				auto &data_all_B = data_all_A[iat2];
				for( auto &data_local_C : data_local_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_local_C.first;
					auto &data_all_C = data_all_B[box2];
					if( data_all_C.c )
						data_all_C += data_local_C.second;
					else
						data_all_C = std::move(data_local_C.second);
				}
			}
		}
	}
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::insert_data(
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> &data_rank,
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all)
{
	for( int is=0; is!=NSPIN; ++is )
	{
		auto &data_rank_is = data_rank[is];
		auto &data_all_is = data_all[is];
		for( auto &data_rank_A : data_rank_is )
		{
			const size_t iat1 = data_rank_A.first;
			auto &data_all_A = data_all_is[iat1];
			for( auto &data_rank_B : data_rank_A.second )
			{
				const size_t iat2 = data_rank_B.first;
				auto &data_all_B = data_all_A[iat2];
				for( auto &data_rank_C : data_rank_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_rank_C.first;
					auto &data_all_C = data_all_B[box2];
					if(data_all_C.c)
						data_all_C += data_rank_C.second.to_matrix();
					else
						data_all_C = data_rank_C.second.to_matrix();
				}
			}
		}
	}
}
