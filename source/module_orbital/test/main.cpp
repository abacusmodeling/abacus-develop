#include "gtest/gtest.h"
int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc,&argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
}
