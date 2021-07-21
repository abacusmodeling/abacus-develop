#include "exx_abfs-parallel-communicate-function.h"
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"

vector<pair<vector<bool>,vector<bool>>>
Exx_Abfs::Parallel::Communicate::Function::get_atom_in_2D_list(const MPI_Comm &mpi_comm)
{
	TITLE("Exx_Abfs::Parallel::Communicate::Functions::get_atom_in_2D_list");
	int comm_sz;	MPI_Comm_size(mpi_comm, &comm_sz);
	
	bool atom_in_2D[ucell.nat*2];
	for(int i=0; i<ucell.nat*2; ++i)
		atom_in_2D[i]=false;
	for(int iwt=0; iwt<GlobalV::NLOCAL; ++iwt)
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