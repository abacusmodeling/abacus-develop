#include "../pw_basis.h"
#include "test_tool.h"

int main(int argc,char **argv)
{

    ModuleBase::Matrix3 latvec = ModuleBase::Matrix3(1,0,0,0,1,0,0,0,1);
    // bool gamma_only = false;
    // double ecut = 10;
    // int nproc, nrank;
    int distribution_type = 1;

    // MPI_Comm_size(POOL_WORLD, &nproc);
    // MPI_Comm_rank(POOL_WORLD, &nrank);

    // ModulePW::PW_Basis pw;
    // pw.initgrids(latvec, ecut);
    // pw.initparameters(gamma_only, ecut, nproc, nrank, distribution_type);
    // pw.distribute_g();

    // for (int ip = 0; ip < nproc; ip++)
    // {
    //     if (nrank == ip)
    //     {
    //         cout<<"ip:  "<<ip<<'\n';
    //         if (ip == 0)
    //         {
    //             cout << "\nnstnz_per  ";
    //             for (int ip = 0; ip < pw.nstot) cout << pw.nstnz_per[ip];
    //             cout << "\nstartnsz_per ";
    //             for (int ip = 0; ip < pw.nstot) cout << pw.nstnz_per[ip];
    //         }
            
    //     }
    // }
    int nproc, myrank;
    int nproc_in_pool, npool, mypool, rank_in_pool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, npool, mypool, rank_in_pool);
    ModulePW::PW_Basis pwtest;


    return 0;
}