//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis.h"
#include "../../module_base/matrix3.h"
#include "test_tool.h"
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
#include <iostream>
#include <complex>

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
    bool gamma_only;
    //--------------------------------------------------
    lat0 = 5;
    ModuleBase::Matrix3 la(1, 0, 0, 0, 1, 0, 0, 0, 1);
    latvec = la;
    wfcecut = 5;
    npool = 1;
    gamma_only = false;
    //--------------------------------------------------
    
    //setup mpi
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, npool, mypool, rank_in_pool);
    
    //init
    pwtest.initgrids(lat0,latvec,4*wfcecut);
    pwtest.initparameters(gamma_only,wfcecut,nproc_in_pool,rank_in_pool,1);
    pwtest.setuptransform();

    int npw = pwtest.npw;
    int nrxx = pwtest.nrxx;
    nx = pwtest.nx;
    ny = pwtest.ny;
    nz = pwtest.nz;
    int nxyz = nx * ny * nz;
    cout<<"FFT: "<<nx<<" "<<ny<<" "<<nz<<endl;
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    
    complex<double> *tmp = NULL;
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
                    double vy = iy -  int(ny/2)+1;
                    double vz = iz -  int(nz/2)+1;
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
        cout << "old method\n";
        ModuleBase::Vector3<double> delta_g(double((int(nx/2)+1))/nx, double((int(ny/2)+1))/ny, double((int(ny/2)+1))/nz); 
        for(int i = 0 ; i < nx*ny ; i+=1)
        {
            int ix = i / (ny * nz);
            int iy = (i - ix*ny*nz)/ nz;
            int iz = i % nz;
            ModuleBase::Vector3<double> real_r(ix, iy, iz);
            double phase_im = delta_g * real_r;
            complex<double> phase(0,ModuleBase::TWO_PI * phase_im);
            tmp[i] /= nxyz;
            tmp[i] *= exp(phase);
            cout<<tmp[i]<<" ";
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
    cout << "new method\n";
    for(int i = 0 ; i < nx*ny ; i+=1)
    {
        cout<<rhor[i]<<" ";
    }
    cout<<endl;

    return 0;
}