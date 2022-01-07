#ifdef __MPI
#include "test_tool.h"
#include "mpi.h"
#endif
#include "fftw3.h"
#include "utest.h"
using namespace std;

int nproc_in_pool, rank_in_pool;

class TestEnv : public testing::Environment 
{
public:
    virtual void SetUp()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"[ SET UP TESTS ]"<<"\033[0m"<<endl;
            cout<<"\033[32m[ "<<nproc_in_pool<<" processors are used."<<" ]\033[0m"<<endl;
        }
    }
    virtual void TearDown()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"[ TEAR DOWN TESTS ]"<<"\033[0m"<<endl;
        }
    }
};



int main(int argc, char **argv) 
{
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
    testing::AddGlobalTestEnvironment(new TestEnv);
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif  
    return result;
}
