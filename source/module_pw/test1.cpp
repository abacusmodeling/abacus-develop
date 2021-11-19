#include "pw_basis.h"

#include "iostream"

int main()
{
    ModuleBase::Matrix3 latvec = ModuleBase::Matrix3(1,0,0,0,1,0,0,0,1);
    bool gamma_only = false;
    double ecut = 10;
    int nproc, nrank;
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
}