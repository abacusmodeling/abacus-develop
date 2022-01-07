#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "../../src_parallel/parallel_global.h"
#include "mpi.h"
#endif
#include "../../module_base/global_function.h"
#include "gtest/gtest.h"
extern int nproc_in_pool,rank_in_pool;

TEST(PWTEST,test1_1)
{
    //--------------------------------------------------
    ModuleBase::Matrix3 latvec(1,0,0,0,1,0,0,0,1);
    bool gamma_only = true;
    double ecut = 100;
    double lat0 = 1;
    int distribution_type = 1;
    //--------------------------------------------------

    ModulePW::PW_Basis pwtest;

    pwtest.initgrids(lat0, latvec, ecut);
    pwtest.initparameters(gamma_only, ecut, nproc_in_pool, rank_in_pool, distribution_type);
    pwtest.distribute_r();
    pwtest.distribute_g();

    int tot_npw = 0;
    int nxy = pwtest.nx * pwtest.ny;
#ifdef __MPI
    MPI_Reduce(&pwtest.npw, &tot_npw, 1, MPI_INT, MPI_SUM, 0, POOL_WORLD);
#else
    tot_npw = pwtest.npw;
#endif

    for (int ip = 0; ip < nproc; ip++)
    {
#ifdef __MPI
        MPI_Barrier(POOL_WORLD);
#endif
        if (rank_in_pool == ip)
        {
            std::cout<<"ip:  "<<ip<<'\n';
            if (ip == 0)
            {
                std::cout<<"nxnynz"<<pwtest.nx<<pwtest.ny<<pwtest.nz<<"\n";
                std::cout << "\nnstnz_per  ";
                for (int ip = 0; ip < nproc; ip++) std::cout << pwtest.nstnz_per[ip] << std::setw(4);
                std::cout << "\nstartnsz_per ";
                for (int ip = 0; ip < nproc; ip++) std::cout << pwtest.startnsz_per[ip] << std::setw(4);
                std::cout << "\nistot2ixy ";
                for (int is = 0; is < pwtest.nstot; is++) std::cout << pwtest.istot2bigixy[is] << std::setw(4);
                std::cout << "\nixy2istot ";
                for (int ixy = 0; ixy < nxy; ixy++) std::cout << pwtest.ixy2istot[ixy] << std::setw(4);
                std::cout << "\nixy2ip ";
                for (int ixy = 0; ixy < nxy; ixy++) std::cout << pwtest.ixy2ip[ixy] << std::setw(4);
                std::cout << "\n";
            }
            std::cout << "ig2isz    ";
            for (int ig = 0; ig < pwtest.npw; ++ig) std::cout << pwtest.ig2isz[ig] << std::setw(4);
            std::cout << "\nis2ixy    ";
            for (int is = 0; is < pwtest.nst; ++is) std::cout << pwtest.is2ixy[is] << std::setw(4);
            std::cout << "\n";
        }
    }
#ifdef __MPI
    MPI_Barrier(POOL_WORLD);
#endif
    pwtest.collect_local_pw();
    for (int ip = 0; ip < nproc; ip++)
    {
#ifdef __MPI
        MPI_Barrier(POOL_WORLD);
#endif
        if (rank_in_pool == ip)
        {
            std::cout<<"ip:  "<<ip<<'\n';
            std::cout << "gg    gdirect    gcar\n";
            for (int ig = 0; ig < pwtest.npw; ++ig)
            {
                std::cout << pwtest.gg[ig] << std::setw(4) << pwtest.gdirect[ig]<< std::setw(4)<< pwtest.gcar[ig] << "\n";
            }
            std::cout << "ig2isz    ";
            for (int ig = 0; ig < pwtest.npw; ++ig) std::cout << pwtest.ig2isz[ig] << std::setw(4);
            std::cout << "\nis2ixy    ";
            for (int is = 0; is < pwtest.nst; ++is) std::cout << pwtest.is2ixy[is] << std::setw(4);
            std::cout << "\n";
        }
    }
}