#include "../pw_basis.h"
#include "test_tool.h"

using namespace std;
int main(int argc,char **argv)
{
    int nproc, myrank;
    int nproc_in_pool, npool, mypool, rank_in_pool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, npool, mypool, rank_in_pool);
    ModulePW::PW_Basis pwtest;
    ModuleBase::Matrix3 G,latvec;
    int nx,ny,nz;  //f*G
    double ggecut;
    double lat0;
    pwtest.initgrids(lat0,latvec,ggecut);
    nx = pwtest.nx;
    ny = pwtest.ny;
    nz = pwtest.nz;
    pwtest.initparameters(false,ggecut,nproc_in_pool,rank_in_pool,1);
    pwtest.distribute_r();
    pwtest.distribute_g();
    int npw = pwtest.npw;
    
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
                    ModuleBase::Vector3<double> v(double(ix),double(iy),double(iz));
                    double modulus = v * (this->GGT * v);
                    if (modulus <= this->ggecut)
                    {
                        tmp[ix*ny*nz + iy*nz + iz]=1.0;
                    }
                }
            }
        }
        fftw_plan pp = fftw_plan_dft_3d(nz,ny,nx,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(pp);
        fftw_free(pp);
    }


    return 0;
}