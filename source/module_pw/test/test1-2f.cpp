//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "../../src_parallel/parallel_global.h"
#include "mpi.h"
#endif
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
#include "pw_test.h"

using namespace std;
TEST_F(PWTEST,test1_2f)
{
    cout<<"dividemthd 1, gamma_only: off, float precision, check fft"<<endl;
    ModulePW::PW_Basis pwtest;
    ModuleBase::Matrix3 latvec;
    int fftnx,fftny,fftnz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    //--------------------------------------------------
    lat0 = 2;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 2, 0, 0, 0, 2);
    latvec = la;
    wfcecut = 10;
    gamma_only = false;
    int distribution_type = 1;
    //--------------------------------------------------
    
    //init
    pwtest.initgrids(lat0,latvec,wfcecut, nproc_in_pool, rank_in_pool);
    //pwtest.initgrids(lat0,latvec,5,7,7);
    pwtest.initparameters(gamma_only,wfcecut,distribution_type);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

    const int npw = pwtest.npw;
    const int nrxx = pwtest.nrxx;
    const int nmaxgr = pwtest.nmaxgr;
    fftnx = pwtest.fftnx;
    fftny = pwtest.fftny;
    fftnz = pwtest.fftnz;
    int nplane = pwtest.nplane;

    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    complex<double> *tmp = new complex<double> [fftnx*fftny*fftnz];
    if(rank_in_pool == 0)
    {
        for(int ix = 0 ; ix < fftnx ; ++ix)
        {
            for(int iy = 0 ; iy < fftny ; ++iy)
            {
                for(int iz = 0 ; iz < fftnz ; ++iz)
                {
                    tmp[ix*fftny*fftnz + iy*fftnz + iz]=0.0;
                    double vx = ix -  int(fftnx/2);
                    double vy = iy -  int(fftny/2);
                    double vz = iz -  int(fftnz/2);
                    ModuleBase::Vector3<double> v(vx,vy,vz);
                    double modulus = v * (GGT * v);
                    if (modulus <= ggecut)
                    {
                        tmp[ix*fftny*fftnz + iy*fftnz + iz]=1.0/(modulus+1) + ModuleBase::IMAG_UNIT / (abs(v.x+1) + 1);
                    }
                }
            }   
        }
        fftw_plan pp = fftw_plan_dft_3d(fftnx,fftny,fftnz,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pp);    
        fftw_destroy_plan(pp); 
        
        ModuleBase::Vector3<double> delta_g(double(int(fftnx/2))/fftnx, double(int(fftny/2))/fftny, double(int(fftny/2))/fftnz); 
        for(int ixy = 0 ; ixy < fftnx * fftny ; ++ixy)
        {
            for(int iz = 0 ; iz < fftnz ; ++iz)
            {
                int ix = ixy / fftny;
                int iy = ixy % fftny;
                ModuleBase::Vector3<double> real_r(ix, iy, iz);
                double phase_im = -delta_g * real_r;
                complex<double> phase(0,ModuleBase::TWO_PI * phase_im);
                tmp[ixy * fftnz + iz] *= exp(phase);
            }
        }
    }
#ifdef __MPI
    MPI_Bcast(tmp,2*fftnx*fftny*fftnz,MPI_DOUBLE,0,POOL_WORLD);
#endif
    
    complex<float> * rhog = new complex<float> [npw];
    complex<float> * rhogr = new complex<float> [nmaxgr];
    complex<float> * rhogout = new complex<float> [npw];
    for(int ig = 0 ; ig < npw ; ++ig)
    {
        rhog[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (abs(pwtest.gdirect[ig].x+1) + 1);
        rhogr[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (abs(pwtest.gdirect[ig].x+1) + 1);
    }    
    complex<float> * rhor = new complex<float> [nrxx];

    pwtest.recip2real(rhog,rhor);//check out-of-place transform

    pwtest.recip2real(rhogr,rhogr);//check in-place transform

    int startiz = pwtest.startz[rank_in_pool];
    for(int ixy = 0 ; ixy < fftnx * fftny ; ++ixy)
    {
        for(int iz = 0 ; iz < nplane ; ++iz)
        {
            EXPECT_NEAR(tmp[ixy * fftnz + startiz + iz].real(),rhor[ixy*nplane+iz].real(),1e-6);
            EXPECT_NEAR(tmp[ixy * fftnz + startiz + iz].imag(),rhor[ixy*nplane+iz].imag(),1e-6);
            EXPECT_NEAR(tmp[ixy * fftnz + startiz + iz].real(),rhogr[ixy*nplane+iz].real(),1e-6);
            EXPECT_NEAR(tmp[ixy * fftnz + startiz + iz].imag(),rhogr[ixy*nplane+iz].imag(),1e-6);
        }
    }

    
    pwtest.real2recip(rhor,rhogout);//check out-of-place transform

    pwtest.real2recip(rhogr,rhogr);//check in-place transform

    for(int ig = 0 ; ig < npw ; ++ig)
    {
        EXPECT_NEAR(rhog[ig].real(),rhogout[ig].real(),1e-6);
        EXPECT_NEAR(rhog[ig].imag(),rhogout[ig].imag(),1e-6);
        EXPECT_NEAR(rhogr[ig].real(),rhogout[ig].real(),1e-6);
        EXPECT_NEAR(rhogr[ig].imag(),rhogout[ig].imag(),1e-6);
    }
  
    delete [] rhog;
    delete [] rhogout;
    delete [] rhor;
    delete []tmp; 
    delete [] rhogr;

    fftw_cleanup();
#ifdef __MIX_PRECISION
    fftwf_cleanup();
#endif
}