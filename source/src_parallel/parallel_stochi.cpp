#include "parallel_stochi.h"
#include "parallel_global.h"
#include "../module_base/global_variable.h"
#include "assert.h"


#ifdef __MPI
void Parallel_Stochi::divide_stogroups(void)
{
    int i=0;
    int j=0;
    if (GlobalV::NPROC < GlobalV::KPAR)
    {
        std::cout<<"\n NPROC=" << GlobalV::NPROC << " NSTOGROUP=" << GlobalV::NSTOGROUP;
        std::cout<<"Error : Too many stogroups !"<<std::endl;
        exit(0);
    }

    // (1) per process in each stogroup
    assert(GlobalV::NPROC%GlobalV::NSTOGROUP==0);
    GlobalV::NPROC_IN_STOGROUP = GlobalV::NPROC/GlobalV::NSTOGROUP;
    GlobalV::MY_STOGROUP = int(GlobalV::MY_RANK / GlobalV::NPROC_IN_STOGROUP);
    GlobalV::RANK_IN_STOGROUP = GlobalV::MY_RANK%GlobalV::NPROC_IN_STOGROUP;

    // (2) per process in each pool
    GlobalV::NPROC_IN_POOL = GlobalV::NPROC_IN_STOGROUP/GlobalV::KPAR;
    if (GlobalV::RANK_IN_STOGROUP < (GlobalV::NPROC_IN_STOGROUP%GlobalV::KPAR)*(GlobalV::NPROC_IN_POOL+1))
    {
        GlobalV::NPROC_IN_POOL++;
        GlobalV::MY_POOL = int(GlobalV::RANK_IN_STOGROUP / GlobalV::NPROC_IN_POOL);
        GlobalV::RANK_IN_POOL = GlobalV::RANK_IN_STOGROUP%GlobalV::NPROC_IN_POOL;
    }
    else
    {
        GlobalV::MY_POOL = int( (GlobalV::RANK_IN_STOGROUP-GlobalV::NPROC_IN_STOGROUP%GlobalV::KPAR) / GlobalV::NPROC_IN_POOL);
        GlobalV::RANK_IN_POOL = (GlobalV::RANK_IN_STOGROUP-GlobalV::NPROC_IN_STOGROUP%GlobalV::KPAR)%GlobalV::NPROC_IN_POOL;
    }
    



    int key = 1;
    MPI_Comm_split(MPI_COMM_WORLD,GlobalV::MY_STOGROUP,key,&STO_WORLD);

    //========================================================
    // MPI_Comm_Split: Creates new communicators based on
    // colors(2nd parameter) and keys(3rd parameter)
    // Note: The color must be non-negative or MPI_UNDEFINED.
    //========================================================
    MPI_Comm_split(STO_WORLD,GlobalV::MY_POOL,key,&POOL_WORLD);
	int color = GlobalV::MY_RANK % GlobalV::NPROC_IN_STOGROUP;
	MPI_Comm_split(MPI_COMM_WORLD, color, key, &PARAPW_WORLD);

    return;
}
#endif

