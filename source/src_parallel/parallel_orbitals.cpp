#include "parallel_orbitals.h"

Parallel_Orbitals::Parallel_Orbitals()
{
    out_hs = 0;

    trace_loc_row = new int[1];
    trace_loc_col = new int[1];

    sender_index_size=1;
	sender_local_index = new int[1];
    sender_size_process = new int[1];
    sender_displacement_process = new int[1];
    sender_size=1;
    sender_buffer=new double[1];

    receiver_index_size=1;
    receiver_global_index = new int[1];
    receiver_size_process = new int[1];
    receiver_displacement_process = new int[1];
    receiver_size=1;
    receiver_buffer=new double[1];
}

Parallel_Orbitals::~Parallel_Orbitals()
{
    delete[] trace_loc_row;
    delete[] trace_loc_col;
    delete[] sender_local_index;
    delete[] sender_size_process;
    delete[] sender_displacement_process;
    delete[] sender_buffer;

    delete[] receiver_global_index;
    delete[] receiver_size_process;
    delete[] receiver_displacement_process;
    delete[] receiver_buffer;
}

bool Parallel_Orbitals::in_this_processor(const int &iw1_all, const int &iw2_all)
{
    if (trace_loc_row[iw1_all] == -1) return false;
    else if (trace_loc_col[iw2_all] == -1) return false;
    return true;
}

void Parallel_Orbitals::set_trace(void)
{
    TITLE("Parallel_Orbitals","set_trace");
    assert(GlobalV::NLOCAL>0);

    delete[] trace_loc_row;
    delete[] trace_loc_col;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"trace_loc_row dimension",GlobalV::NLOCAL);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"trace_loc_col dimension",GlobalV::NLOCAL);

    trace_loc_row = new int[GlobalV::NLOCAL];
    trace_loc_col = new int[GlobalV::NLOCAL];
    // mohan update 2011-04-07
    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        trace_loc_row[i] = -1;
        trace_loc_col[i] = -1;
    }

    ModuleBase::Memory::record("Parallel_Orbitals","trace_loc_row",GlobalV::NLOCAL,"int");
    ModuleBase::Memory::record("Parallel_Orbitals","trace_loc_col",GlobalV::NLOCAL,"int");

    if(GlobalV::KS_SOLVER=="lapack"
    || GlobalV::KS_SOLVER=="cg"
    || GlobalV::KS_SOLVER=="dav") //xiaohui add 2013-09-02
	{
		std::cout << " common settings for trace_loc_row and dtraace_loc_col " << std::endl;
		for (int i=0; i<GlobalV::NLOCAL; i++)
		{
			trace_loc_row[i] = i;
			trace_loc_col[i] = i;
		}
		this->nrow = GlobalV::NLOCAL;
		this->ncol = GlobalV::NLOCAL;
	}
#ifdef __MPI
    else if(GlobalV::KS_SOLVER=="scalpack" || GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="hpseps" 
		|| GlobalV::KS_SOLVER=="selinv" || GlobalV::KS_SOLVER=="scalapack_gvx") //xiaohui add 2013-09-02
    {
        // GlobalV::ofs_running << " nrow=" << nrow << std::endl;
        for (int irow=0; irow< this->nrow; irow++)
        {
            int global_row = MatrixInfo.row_set[irow];
            trace_loc_row[global_row] = irow;
			// GlobalV::ofs_running << " global_row=" << global_row 
			// << " trace_loc_row=" << trace_loc_row[global_row] << std::endl;
        }

        // GlobalV::ofs_running << " ncol=" << ncol << std::endl;
        for (int icol=0; icol< this->ncol; icol++)
        {
            int global_col = MatrixInfo.col_set[icol];
            trace_loc_col[global_col] = icol;
			// GlobalV::ofs_running << " global_col=" << global_col 
			// << " trace_loc_col=" << trace_loc_col[global_col] << std::endl;
        }
    }
#endif
    else 
    {
        std::cout << " Parallel Orbial, GlobalV::DIAGO_TYPE = " << GlobalV::KS_SOLVER << std::endl;
        WARNING_QUIT("Parallel_Orbitals::set_trace","Check GlobalV::KS_SOLVER.");
    }

    //---------------------------
    // print the trace for test.
    //---------------------------
    /*
    GlobalV::ofs_running << " " << std::setw(10) << "GlobalRow" << std::setw(10) << "LocalRow" << std::endl;
    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        GlobalV::ofs_running << " " << std::setw(10) << i << std::setw(10) << trace_loc_row[i] << std::endl;

    }

    GlobalV::ofs_running << " " << std::setw(10) << "GlobalCol" << std::setw(10) << "LocalCol" << std::endl;
    for(int j=0; j<GlobalV::NLOCAL; j++)
    {
        GlobalV::ofs_running << " " << std::setw(10) << j << std::setw(10) << trace_loc_col[j] << std::endl;
    }
    */

    return;
}


