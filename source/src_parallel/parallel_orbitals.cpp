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
    assert(NLOCAL>0);

    delete[] trace_loc_row;
    delete[] trace_loc_col;

    OUT(ofs_running,"trace_loc_row dimension",NLOCAL);
    OUT(ofs_running,"trace_loc_col dimension",NLOCAL);

    trace_loc_row = new int[NLOCAL];
    trace_loc_col = new int[NLOCAL];
    // mohan update 2011-04-07
    for(int i=0; i<NLOCAL; i++)
    {
        trace_loc_row[i] = -1;
        trace_loc_col[i] = -1;
    }

    Memory::record("Parallel_Orbitals","trace_loc_row",NLOCAL,"int");
    Memory::record("Parallel_Orbitals","trace_loc_col",NLOCAL,"int");

    if(KS_SOLVER=="lapack"
    || KS_SOLVER=="cg"
    || KS_SOLVER=="dav") //xiaohui add 2013-09-02
	{
		cout << " common settings for trace_loc_row and dtraace_loc_col " << endl;
		for (int i=0; i<NLOCAL; i++)
		{
			trace_loc_row[i] = i;
			trace_loc_col[i] = i;
		}
		this->nrow = NLOCAL;
		this->ncol = NLOCAL;
	}
#ifdef __MPI
    else if(KS_SOLVER=="scalpack" || KS_SOLVER=="genelpa" || KS_SOLVER=="hpseps" 
		|| KS_SOLVER=="selinv" || KS_SOLVER=="scalapack_gvx") //xiaohui add 2013-09-02
    {
        // ofs_running << " nrow=" << nrow << endl;
        for (int irow=0; irow< this->nrow; irow++)
        {
            int global_row = MatrixInfo.row_set[irow];
            trace_loc_row[global_row] = irow;
			// ofs_running << " global_row=" << global_row 
			// << " trace_loc_row=" << trace_loc_row[global_row] << endl;
        }

        // ofs_running << " ncol=" << ncol << endl;
        for (int icol=0; icol< this->ncol; icol++)
        {
            int global_col = MatrixInfo.col_set[icol];
            trace_loc_col[global_col] = icol;
			// ofs_running << " global_col=" << global_col 
			// << " trace_loc_col=" << trace_loc_col[global_col] << endl;
        }
    }
#endif
    else 
    {
        cout << " Parallel Orbial, DIAGO_TYPE = " << KS_SOLVER << endl;
        WARNING_QUIT("Parallel_Orbitals::set_trace","Check KS_SOLVER.");
    }

    //---------------------------
    // print the trace for test.
    //---------------------------
    /*
    ofs_running << " " << setw(10) << "GlobalRow" << setw(10) << "LocalRow" << endl;
    for(int i=0; i<NLOCAL; i++)
    {
        ofs_running << " " << setw(10) << i << setw(10) << trace_loc_row[i] << endl;

    }

    ofs_running << " " << setw(10) << "GlobalCol" << setw(10) << "LocalCol" << endl;
    for(int j=0; j<NLOCAL; j++)
    {
        ofs_running << " " << setw(10) << j << setw(10) << trace_loc_col[j] << endl;
    }
    */

    return;
}


