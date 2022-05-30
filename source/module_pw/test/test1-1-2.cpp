#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "../../src_parallel/parallel_global.h"
#include "mpi.h"
#endif
#include "../../module_base/global_function.h"
#include "../../module_base/constants.h"
#include "pw_test.h"
extern int nproc_in_pool,rank_in_pool;
using namespace std;

TEST_F(PWTEST,test1_1_2)
{
    cout<<"dividemthd 1, gamma_only: on, check gcar,gdirect,gg,istot2bigixy,ig2isz"<<endl;
    //--------------------------------------------------
    ModuleBase::Matrix3 latvec(1,2,0,2,1,1,0,0,5);
    bool gamma_only = true;
    double wfcecut = 70;
    double lat0 = 4;
    int distribution_type = 1;
    //--------------------------------------------------

    ModulePW::PW_Basis pwtest;

    pwtest.initgrids(lat0, latvec, wfcecut,nproc_in_pool, rank_in_pool);
    pwtest.initparameters(gamma_only, wfcecut, distribution_type);
    pwtest.setuptransform();
    pwtest.collect_local_pw();
    pwtest.collect_uniqgg();
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;

    //ref
    const int totnpw_ref = 5012;
    const int totnst_ref = 176;
    const int nx_ref = 24;
    const int bigny_ref = 27;
    const int ny_ref = 14;
    const int nz_ref = 54;
    //some results for different number of processors
    int npw_per_ref[12][12]={
        {5012},
        {2506,2506},
        {1672,1670,1670},
        {1253,1253,1251,1255},
        {1002,1002,1002,1003,1003},
        {835,837,835,836,835,834},
        {716,717,716,716,716,716,715},
        {627,627,628,626,626,627,625,626},
        {556,559,555,555,555,557,559,559,557},
        {502,502,501,499,498,502,502,503,503,500},
        {457,457,454,458,458,454,456,455,452,455,456},
        {419,419,419,418,417,418,416,416,416,420,416,418}
    };
    int nst_per_ref[12][12]={
        {176},
        {88,88},
        {59,59,58},
        {44,44,44,44},
        {35,35,35,36,35},
        {29,30,30,29,29,29},
        {25,25,25,25,26,25,25},
        {22,22,22,22,22,22,22,22},
        {20,19,19,19,19,20,20,20,20},
        {17,17,18,18,17,17,18,18,18,18},
        {16,16,16,16,16,16,16,16,16,16,16},
        {15,15,15,15,15,15,14,14,14,15,14,15}
    };
    int *npw_per;
    if(rank_in_pool == 0)
    {
        npw_per = new int [nproc_in_pool];
    }
#ifdef __MPI
    MPI_Gather(&pwtest.npw,1,MPI_INT,npw_per,1,MPI_INT,0,POOL_WORLD);
#else
    if(rank_in_pool == 0) npw_per[0] = pwtest.npw;
#endif
    if(rank_in_pool == 0)
    {
        if(nproc_in_pool <= 12)
        {
            for(int ip = 0 ; ip < nproc_in_pool ; ++ip)
            {
                ASSERT_EQ(npw_per_ref[nproc_in_pool-1][ip], npw_per[ip]);
                ASSERT_EQ(nst_per_ref[nproc_in_pool-1][ip], pwtest.nst_per[ip]);
            }
        }
        else
        {
            cout<<"Please use mpi processors no more than 12."<<endl;
        }
        delete []npw_per;
    }
    

    //results
    int tot_npw = 0;
#ifdef __MPI
    MPI_Allreduce(&pwtest.npw, &tot_npw, 1, MPI_INT, MPI_SUM, POOL_WORLD);
#else
    tot_npw = pwtest.npw;
#endif
    ASSERT_EQ(pwtest.nx, nx_ref);
    ASSERT_EQ(pwtest.ny, ny_ref);
    ASSERT_EQ(pwtest.bigny, bigny_ref);
    ASSERT_EQ(pwtest.nz, nz_ref);
    ASSERT_EQ(tot_npw, totnpw_ref);
    ASSERT_EQ(pwtest.npwtot, totnpw_ref);
    ASSERT_EQ(pwtest.nstot,totnst_ref);
    ASSERT_EQ(pwtest.bignxyz, nx_ref*bigny_ref*nz_ref);


    int *tmpx = new int[pwtest.nx*pwtest.ny*pwtest.nz];
    int *tmpy = new int[pwtest.nx*pwtest.ny*pwtest.nz];
    int *tmpz = new int[pwtest.nx*pwtest.ny*pwtest.nz];
    ModuleBase::GlobalFunc::ZEROS(tmpx,pwtest.nx*pwtest.ny*pwtest.nz);
    ModuleBase::GlobalFunc::ZEROS(tmpy,pwtest.nx*pwtest.ny*pwtest.nz);
    ModuleBase::GlobalFunc::ZEROS(tmpz,pwtest.nx*pwtest.ny*pwtest.nz);
    
    int * startnst = new int [nproc_in_pool];
    startnst[0] = 0;
    for(int ip = 1 ; ip < nproc_in_pool; ++ip)
    {
        startnst[ip] = startnst[ip-1] + pwtest.nst_per[ip-1];
    }

    for(int ig = 0 ; ig < pwtest.npw; ++ig)
    {
        int istot = pwtest.ig2isz[ig] / pwtest.nz + startnst[rank_in_pool];
        // int is = pwtest.ig2isz[ig] / pwtest.nz;
        int iz = pwtest.ig2isz[ig] % pwtest.nz;
        int iy = pwtest.istot2bigixy[istot] % pwtest.bigny;
        int ix = pwtest.istot2bigixy[istot] / pwtest.bigny;
        // int iy = pwtest.is2ixy[is] % pwtest.ny;
        // int ix = pwtest.is2ixy[is] / pwtest.ny;

        tmpx[iz+(iy+ix*pwtest.ny)*pwtest.nz] = int(pwtest.gdirect[ig].x);
        tmpy[iz+(iy+ix*pwtest.ny)*pwtest.nz] = int(pwtest.gdirect[ig].y);
        tmpz[iz+(iy+ix*pwtest.ny)*pwtest.nz] = int(pwtest.gdirect[ig].z);
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,tmpx,pwtest.nxyz,MPI_INT,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,tmpy,pwtest.nxyz,MPI_INT,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,tmpz,pwtest.nxyz,MPI_INT,MPI_SUM,POOL_WORLD);
#endif
    if(rank_in_pool==0)
    {
        for(int iz = 0 ; iz < pwtest.nz; ++iz)
        {
            for(int iy = 0 ; iy < pwtest.ny ; ++iy)
            {
                for(int ix = 0 ; ix < pwtest.nx ; ++ix)
                {
                    ModuleBase::Vector3<double> f;
                    f.x = ix;
                    f.y = iy;
                    f.z = iz;
                    if(iz >= int(pwtest.nz/2) +1) f.z -= pwtest.nz;
                    if(ix >= int(pwtest.nx/2) +1) f.x -= pwtest.nx;
                    double modulus = f * (GGT * f);
                    if (modulus <= ggecut)
                    {
                        EXPECT_EQ(tmpx[iz + iy*pwtest.nz + ix*pwtest.ny*pwtest.nz], int(f.x));
                        EXPECT_EQ(tmpy[iz + iy*pwtest.nz + ix*pwtest.ny*pwtest.nz], int(f.y));
                        EXPECT_EQ(tmpz[iz + iy*pwtest.nz + ix*pwtest.ny*pwtest.nz], int(f.z));
                    }
                    
                }
            }
        }
    }
    for(int ig = 0 ;ig < pwtest.npw ; ++ig)
    {
        ModuleBase::Vector3<double> f;
        f.x = pwtest.gdirect[ig].x;
        f.y = pwtest.gdirect[ig].y;
        f.z = pwtest.gdirect[ig].z;
        ModuleBase::Vector3<double> gcar;
        gcar = f * G;
        double modulus = f*GGT*f;
        EXPECT_NEAR(gcar.x,pwtest.gcar[ig].x,1e-6);
        EXPECT_NEAR(gcar.y,pwtest.gcar[ig].y,1e-6);
        EXPECT_NEAR(gcar.z,pwtest.gcar[ig].z,1e-6);
        EXPECT_NEAR(modulus,pwtest.gg[ig],1e-6);
        EXPECT_NEAR(pwtest.gg[ig], pwtest.gg_uniq[pwtest.ig2igg[ig]],1e-8);
    }
    for(int igg = 1 ; igg < pwtest.ngg ; ++igg)
    {
        EXPECT_GT(pwtest.gg_uniq[igg], pwtest.gg_uniq[igg-1]);
    }
    if(pwtest.ig_gge0 >= 0) EXPECT_NEAR(0.0, pwtest.gg[pwtest.ig_gge0], 1e-8);
    delete [] startnst;
    delete [] tmpx;
    delete [] tmpy;
    delete [] tmpz;
}