#ifdef __DEEPKS

#include "parallel_deepks.h"

namespace GlobalC
{
Parallel_deepks ParaD;
}

Parallel_deepks::Parallel_deepks()
{
    this->trace_loc_orb=new int[1];
    this->norb_local=new int[1];
}

Parallel_deepks::~Parallel_deepks()
{
    delete[] this->trace_loc_orb;
    delete[] this->norb_local;
}

void Parallel_deepks::set_nlocal(
#ifdef __MPI
		MPI_Comm MPI_WORLD
#endif
)
{
    ModuleBase::TITLE("Parallel_deepks","cal_nlocal");
    assert(GlobalV::NLOCAL>0);

    delete[] this->norb_local;
    
    int size;
#ifdef __MPI
    assert(GlobalV::NPROC>0);
    size = GlobalV::NPROC;
#else
    size = 1;
#endif
    this->norb_local=new int[size];

#ifdef __MPI
    int loc_size_base=GlobalV::NLOCAL/GlobalV::NPROC;
    if(loc_size_base==0)
    {
        GlobalV::ofs_warning << " loc_size=0" << " in proc " << GlobalV::MY_RANK+1 << std::endl;
        ModuleBase::WARNING_QUIT("Parallel_deepks::set_nlocal","NLOCAL < GlobalV::NPROC");
    }

    int sum_nlocal = 0;
    for(int rank=0;rank<GlobalV::NPROC;rank++)
    {
        int loc_size=loc_size_base;
        if (rank<GlobalV::NLOCAL%GlobalV::NPROC) loc_size+=1;
        this->norb_local[rank] = loc_size;
        sum_nlocal += loc_size;
        //GlobalV::ofs_running << " loc_size: " << loc_size << " in proc " << rank << std::endl;
        //GlobalV::ofs_running << " sum_nlocal: " << sum_nlocal << " in proc " << GlobalV::MY_RANK+1 << std::endl;
    }

    //check if sum of norb_local = tot. nlocal
    assert(sum_nlocal == GlobalV::NLOCAL);
#else
    this->norb_local[0] = GlobalV::NLOCAL;
#endif
}

void Parallel_deepks::set_loc_orb
(
#ifdef __MPI
		MPI_Comm MPI_WORLD
#endif
)
{
	ModuleBase::TITLE("Parallel_deepks","set_loc_orb");
	assert(GlobalV::NLOCAL>0);

    delete[] this->trace_loc_orb;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"trace_loc_orb dimension",GlobalV::NLOCAL);

    this->trace_loc_orb = new int[GlobalV::NLOCAL];
    ModuleBase::Memory::record("Parallel_deepks","trace_loc_orb",GlobalV::NLOCAL,"int");

    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        trace_loc_orb[i] = -1;
    }

#ifdef __MPI
    int norb = this->norb_local[GlobalV::MY_RANK];

    int norb_start = 0;
    //this is the global index of first local orb on current processor
    for(int i=0;i<GlobalV::MY_RANK;i++)
    {
        norb_start += this->norb_local[i];
    }
    
    for (int iorb=0; iorb< norb; iorb++)
    {
        int global_orb = iorb + norb_start;
        this->trace_loc_orb[global_orb] = iorb;
    }
    
    //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"norb:",norb);
    //for(int i=0;i<GlobalV::NLOCAL;i++)
    //{
        //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"trace:",this->trace_loc_orb[i]);
    //}
#else
    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        trace_loc_orb[i] = i;
    }
#endif
}

void Parallel_deepks::allsum_deepks(
	int inlmax, //first dimension
	int ndim, //second dimension
	double** mat)
{
    for(int inl=0;inl<inlmax;inl++)
    {
        Parallel_Reduce::reduce_double_all(mat[inl],ndim);
    }
}
#endif