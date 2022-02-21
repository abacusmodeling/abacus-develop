#include "parallel_orbitals.h"
#include "../module_base/memory.h"
#ifdef __MPI
extern "C"
{
    #include "../module_base/blacs_connector.h"
	#include "../module_base/scalapack_connector.h"
}
#endif

Parallel_Orbitals::Parallel_Orbitals()
{
    loc_sizes = new int[1];
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
    ModuleBase::TITLE("Parallel_Orbitals","set_trace");
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
        ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_trace","Check ks_solver.");
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

#ifdef __MPI
inline int cart2blacs(
	MPI_Comm comm_2D,
	int nprows,
	int npcols,
	int N,
	int nblk,
	int lld,
	int *desc)
{
    int my_blacs_ctxt;
    int myprow, mypcol;
    int *usermap=new int[nprows*npcols];
    int info=0;
    for(int i=0; i<nprows; ++i)
    {
        for(int j=0; j<npcols; ++j)
        {
            int pcoord[2]={i, j};
            MPI_Cart_rank(comm_2D, pcoord, &usermap[i+j*nprows]);
        }
    }
    MPI_Fint comm_2D_f = MPI_Comm_c2f(comm_2D);
    Cblacs_get(comm_2D_f, 0, &my_blacs_ctxt);
    Cblacs_gridmap(&my_blacs_ctxt, usermap, nprows, nprows, npcols);
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    delete[] usermap;
    int ISRC=0;
    descinit_(desc, &N, &N, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);

    return my_blacs_ctxt;
}
#endif

void Parallel_Orbitals::divide_HS_2d
(
#ifdef __MPI
	MPI_Comm DIAG_WORLD
#endif
)
{
	ModuleBase::TITLE("Parallel_Orbitals","divide_HS_2d");
	assert(GlobalV::NLOCAL>0);
	assert(GlobalV::DSIZE>0);

#ifdef __MPI
	DIAG_HPSEPS_WORLD=DIAG_WORLD;
#endif

	if(GlobalV::DCOLOR!=0) return; // mohan add 2012-01-13

	// get the 2D index of computer.
	this->dim0 = (int)sqrt((double)GlobalV::DSIZE); //mohan update 2012/01/13
	//while (GlobalV::NPROC_IN_POOL%dim0!=0)
	while (GlobalV::DSIZE%dim0!=0)
	{
		this->dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	this->dim1=GlobalV::DSIZE/dim0;

	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dim0",dim0);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dim1",dim1);

#ifdef __MPI
	// mohan add 2011-04-16
	if(GlobalV::NB2D==0)
	{
		if(GlobalV::NLOCAL>0) this->nb = 1;
		if(GlobalV::NLOCAL>500) this->nb = 32;
		if(GlobalV::NLOCAL>1000) this->nb = 64;
	}
	else if(GlobalV::NB2D>0)
	{
		this->nb = GlobalV::NB2D; // mohan add 2010-06-28
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nb2d",nb);

	this->set_parameters();

	// call mpi_creat_cart
	this->mpi_creat_cart(&this->comm_2D,this->dim0,this->dim1);

	// call mat_2d
	this->mat_2d(this->comm_2D, GlobalV::NLOCAL, GlobalV::NLOCAL, this->nb, this->MatrixInfo);

	// mohan add 2010-06-29
	this->nrow = this->MatrixInfo.row_num;
	this->ncol = this->MatrixInfo.col_num;
	this->nloc = MatrixInfo.col_num * MatrixInfo.row_num;

	// init blacs context for genelpa
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        blacs_ctxt=cart2blacs(comm_2D, dim0, dim1, GlobalV::NLOCAL, nb, nrow, desc);
    }
#else // single processor used.
	this->nb = GlobalV::NLOCAL;
	this->nrow = GlobalV::NLOCAL;
	this->ncol = GlobalV::NLOCAL;
	this->nloc = GlobalV::NLOCAL * GlobalV::NLOCAL;
	this->set_parameters();
	MatrixInfo.row_b = 1;
	MatrixInfo.row_num = GlobalV::NLOCAL;
	delete[] MatrixInfo.row_set;
	MatrixInfo.row_set = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		MatrixInfo.row_set[i]=i;
	}
	MatrixInfo.row_pos=0;

	MatrixInfo.col_b = 1;
	MatrixInfo.col_num = GlobalV::NLOCAL;
	delete[] MatrixInfo.col_set;
	MatrixInfo.col_set = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		MatrixInfo.col_set[i]=i;
	}
	MatrixInfo.col_pos=0;
#endif

	assert(nloc>0);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MatrixInfo.row_num",MatrixInfo.row_num);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MatrixInfo.col_num",MatrixInfo.col_num);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nloc",nloc);
	return;
}
