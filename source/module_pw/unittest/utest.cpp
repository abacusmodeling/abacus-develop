#include "gtest/gtest.h"
#ifdef __MPI
#include "test_tool.h"
#include "mpi.h"
#endif
#include "fftw3.h"
extern int nproc_in_pool, rank_in_pool;
int nproc_in_pool, rank_in_pool;

int main(int argc, char **argv) {
    int nproc, myrank;
    int npool, mypool;
    npool = 1;
#ifdef __MPI
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, npool, mypool, rank_in_pool);
#else
    nproc = nproc_in_pool = npool = 1;
    myrank = mypool = rank_in_pool = 0;
#endif
    int result = 0;
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif  
    return result;
}
