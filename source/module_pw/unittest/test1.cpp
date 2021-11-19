#include "../pw_basis.h"
#include "test_tool.h"

int main(int argc,char **argv)
{
    int nproc, myrank;
    int nproc_in_pool, npool, mypool, rank_in_pool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, npool, mypool, rank_in_pool);
    ModulePW::PW_Basis pwtest;

    return 0;
}