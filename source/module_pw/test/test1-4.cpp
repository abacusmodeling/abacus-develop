//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis_k.h"
#ifdef __MPI
#include "test_tool.h"
#include "../../src_parallel/parallel_global.h"
#include "mpi.h"
#endif
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
#include "pw_test.h"

using namespace std;
TEST_F(PWTEST,test1_4)
{
    cout<<"dividemthd 1, gamma_only: off, double precision, 2 kpoints, check fft"<<endl;
    ModulePW::PW_Basis_K pwtest;
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    ModuleBase::Vector3<double> *kvec_d;
    int nks;
    //--------------------------------------------------
    lat0 = 2;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 2, 0, 0, 0, 2);
    nks = 2;
    kvec_d = new ModuleBase::Vector3<double>[nks];
    kvec_d[0].set(0,0,0.5);
    kvec_d[1].set(0.5,0.5,0.5);
    latvec = la;
    wfcecut = 10;
    gamma_only = false;
    int distribution_type = 1;
    //--------------------------------------------------
    //init
    pwtest.initgrids(lat0,latvec,wfcecut, nproc_in_pool, rank_in_pool);
    //pwtest.initgrids(lat0,latvec,5,7,7);
    pwtest.initparameters(gamma_only,wfcecut,nks,kvec_d,distribution_type);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

    int nrxx = pwtest.nrxx;
    nx = pwtest.nx;
    ny = pwtest.bigny;
    nz = pwtest.nz;
    int nplane = pwtest.nplane;

    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    complex<double> *tmp = new complex<double> [nx*ny*nz];
    complex<double> * rhor = new complex<double> [nrxx];
    for(int ik  = 0; ik < nks; ++ik)
    {
        int npw = pwtest.npwk[ik];
        int npwk_max = pwtest.npwk_max;
        if(rank_in_pool == 0)
        {
            ModuleBase::Vector3<double> kk = kvec_d[ik];
            for(int ix = 0 ; ix < nx ; ++ix)
            {
                for(int iy = 0 ; iy < ny ; ++iy)
                {
                    for(int iz = 0 ; iz < nz ; ++iz)
                    {
                        tmp[ix*ny*nz + iy*nz + iz]=0.0;
                        double vx = ix -  int(nx/2);
                        double vy = iy -  int(ny/2);
                        double vz = iz -  int(nz/2);
                        ModuleBase::Vector3<double> v(vx,vy,vz);
                        double modulus = v * (GGT * v);
                        double modulusgk = (v+kk) * (GGT * (v+kk));
                        if (modulusgk <= ggecut)
                        {
                            tmp[ix*ny*nz + iy*nz + iz]=1.0/(modulus+1) + ModuleBase::IMAG_UNIT / (abs(v.x+1) + 1);
                        }
                    }
                }   
            }
            fftw_plan pp = fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(pp);    
            fftw_destroy_plan(pp); 

            ModuleBase::Vector3<double> delta_g(double(int(nx/2))/nx, double(int(ny/2))/ny, double(int(ny/2))/nz); 
            for(int ixy = 0 ; ixy < nx * ny ; ++ixy)
            {
                for(int iz = 0 ; iz < nz ; ++iz)
                {
                    int ix = ixy / ny;
                    int iy = ixy % ny;
                    ModuleBase::Vector3<double> real_r(ix, iy, iz);
                    double phase_im = -delta_g * real_r;
                    complex<double> phase(0,ModuleBase::TWO_PI * phase_im);
                    tmp[ixy * nz + iz] *= exp(phase);
                }
            }
        }
#ifdef __MPI
        MPI_Bcast(tmp,2*nx*ny*nz,MPI_DOUBLE,0,POOL_WORLD);
#endif
        complex<double> * rhog = new complex<double> [npw];
        complex<double> * rhogout = new complex<double> [npw];
        for(int ig = 0 ; ig < npw ; ++ig)
        {
            rhog[ig] = 1.0/(pwtest.gg[ik*npwk_max + ig]+1) + ModuleBase::IMAG_UNIT / (abs(pwtest.gdirect[ik*npwk_max + ig].x+1) + 1);
        }    

        pwtest.recip2real(rhog,rhor,ik);
        int startiz = pwtest.startz[rank_in_pool];
        for(int ixy = 0 ; ixy < nx * ny ; ++ixy)
        {
            for(int iz = 0 ; iz < nplane ; ++iz)
            {
                EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhor[ixy*nplane+iz].real(),1e-6);
                EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhor[ixy*nplane+iz].imag(),1e-6);
            }
        }

        pwtest.real2recip(rhor,rhogout,ik);
        for(int ig = 0 ; ig < npw ; ++ig)
        {
            EXPECT_NEAR(rhog[ig].real(),rhogout[ig].real(),1e-6);
            EXPECT_NEAR(rhog[ig].imag(),rhogout[ig].imag(),1e-6);
        }


        delete [] rhog;
        delete [] rhogout;
    }
    delete []tmp; 
    delete [] rhor;
    delete[] kvec_d;
    fftw_cleanup();
#ifdef __MIX_PRECISION
    fftwf_cleanup();
#endif
}