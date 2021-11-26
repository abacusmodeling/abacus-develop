#include "../pw_basis.h"
#include "../../module_base/matrix3.h"
#include "test_tool.h"
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
#include <iostream>

using namespace std;
int main(int argc,char **argv)
{
    int nproc, myrank;
    int nproc_in_pool, npool, mypool, rank_in_pool;
    ModulePW::PW_Basis pwtest;
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    //--------------------------------------------------
    lat0 = 5;
    ModuleBase::Matrix3 la(1, 0, 0, 0, 1, 0, 0, 0, 1);
    latvec = la;
    wfcecut = 5;
    npool = 1;
    //--------------------------------------------------
    
    //setup mpi
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, npool, mypool, rank_in_pool);
    
    //init
    pwtest.initgrids(lat0,latvec,4*wfcecut);
    pwtest.initparameters(false,wfcecut,nproc_in_pool,rank_in_pool,1);
    pwtest.distribute();

    int npw = pwtest.npw;
    int nrxx = pwtest.nrxx;
    nx = pwtest.nx;
    ny = pwtest.ny;
    nz = pwtest.nz;
    cout<<"FFT: "<<nx<<" "<<ny<<" "<<nz<<endl;
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    
    if(myrank == 0)
    {
        complex<double> *tmp = new complex<double> [nx*ny*nz];
        for(int ix = 0 ; ix < nx ; ++ix)
        {
            for(int iy = 0 ; iy < ny ; ++iy)
            {
                for(int iz = 0 ; iz < nz ; ++iz)
                {
                    tmp[ix*ny*nz + iy*nz + iz]=0.0;
                    double vx = ix -  int(nx/2)+1;
                    double vy = iy -  int(nx/2)+1;
                    double vz = iz -  int(nx/2)+1;
                    ModuleBase::Vector3<double> v(vx,vy,vz);
                    double modulus = v * (GGT * v);
                    if (modulus <= ggecut)
                    {
                        tmp[ix*ny*nz + iy*nz + iz]=1.0;
                    }
                }
            }
        }
        fftw_plan pp = fftw_plan_dft_3d(nz,ny,nx,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pp);
        fftw_free(pp);
        
        //output
        for(int i = 0 ; i < nx*ny ; i+=4)
        {
            cout<<tmp[nx*ny]<<" ";
        }
        cout<<endl;
    }
    
    complex<double> * rhog = new complex<double> [npw];
    for(int ig = 0 ; ig < npw ; ++ig)
    {
        rhog[ig] = 1.0;
    }    
    complex<double> * rhor = new complex<double> [nrxx];
    pwtest.recip2real(rhog,rhor);
    for(int i = 0 ; i < nx*ny ; i+=4)
    {
        cout<<rhor[nx*ny]<<" ";
    }
    cout<<endl;


    return 0;
}