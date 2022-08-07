#include "pw_basis.h"
#include "../module_base/global_function.h"
#include "../module_base/timer.h"
#include "typeinfo"
namespace ModulePW
{
/// 
/// in: (nplane,fftny,fftnx) out:(nz,nst)
/// in and out should be in different places 
/// in[] will be changed
/// 
template<typename T>
void PW_Basis:: gatherp_scatters(std::complex<T> *in, std::complex<T> *out)
{
    ModuleBase::timer::tick(this->classname, "gatherp_scatters");
    
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
        std::complex<T> * outp, *inp;
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = this->istot2ixy[is];
            //int ixy = (ixy / fftny)*ny + ixy % fftny;
            outp = &out[is*nz];
            inp = &in[ixy*nz];
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        ModuleBase::timer::tick(this->classname, "gatherp_scatters");
        return;
    }
#ifdef __MPI
    //change (nplane fftnxy) to (nplane,nstot)
    // Hence, we can send them at one time.
    std::complex<T> * outp, *inp;  
	for (int istot = 0;istot < nstot; ++istot)
	{
		int ixy = this->istot2ixy[istot];
        //int ixy = (ixy / fftny)*ny + ixy % fftny;
        outp = &out[istot*nplane];
        inp = &in[ixy*nplane];
		for (int iz = 0; iz < nplane; ++iz)
		{
			outp[iz] = inp[iz];
		}
	}

    //exchange data
    //(nplane,nstot) to (numz[ip],ns, poolnproc)
    if(typeid(T) == typeid(double))
	    MPI_Alltoallv(out, numr, startr, mpi_dcomplex, in, numg, startg, mpi_dcomplex, this->pool_world);
    else if(typeid(T) == typeid(float))
        MPI_Alltoallv(out, numr, startr, MPI_COMPLEX, in, numg, startg, MPI_COMPLEX, this->pool_world);
    // change (nz,ns) to (numz[ip],ns, poolnproc)
    std::complex<T> * outp0, *inp0;
    for (int ip = 0; ip < this->poolnproc ;++ip)
	{
        int nzip = this->numz[ip];
        outp0 = &out[startz[ip]];
        inp0 = &in[startg[ip]];
		for (int is = 0; is < this->nst; ++is)
		{
            outp = &outp0[is * nz];
            inp = &inp0[is * nzip ];
			for (int izip = 0; izip < nzip; ++izip)
			{
				outp[izip] = inp[izip];
			}
		}
	}
   
#endif
    ModuleBase::timer::tick(this->classname, "gatherp_scatters");
    return;
}

/// 
/// in: (nz,nst) out:(nplane,fftny,fftnx)
/// in and out should be in different places 
/// in[] will be changed
/// 
template<typename T>
void PW_Basis:: gathers_scatterp(std::complex<T> *in, std::complex<T> *out)
{
    ModuleBase::timer::tick(this->classname, "gathers_scatterp");
    
    if(this->poolnproc == 1) //In this case nrxx=fftnx*fftny*nz, nst = nstot, 
    {
        ModuleBase::GlobalFunc::ZEROS(out, this->nrxx);
        std::complex<T> * outp, *inp;
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = istot2ixy[is];
            //int ixy = (ixy / fftny)*ny + ixy % fftny;
            outp = &out[ixy*nz];
            inp = &in[is*nz];
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        ModuleBase::timer::tick(this->classname, "gathers_scatterp");
        return;
    }
#ifdef __MPI
    // change (nz,ns) to (numz[ip],ns, poolnproc)
    // Hence, we can send them at one time. 
    std::complex<T> * outp, *inp; 
    std::complex<T> * outp0, *inp0;
    for (int ip = 0; ip < this->poolnproc ;++ip)
	{
        int nzip = this->numz[ip];
        outp0 = &out[startg[ip]];
        inp0 = &in[startz[ip]];
		for (int is = 0; is < this->nst; ++is)
		{
            outp = &outp0[is * nzip];
            inp = &inp0[is * nz ];
			for (int izip = 0; izip < nzip; ++izip)
			{
				outp[izip] = inp[izip];
			}
		}
	}

	//exchange data
    //(numz[ip],ns, poolnproc) to (nplane,nstot)
    if(typeid(T) == typeid(double))
	    MPI_Alltoallv(out, numg, startg, mpi_dcomplex, in, numr, startr, mpi_dcomplex, this->pool_world);
    else if(typeid(T) == typeid(float))
        MPI_Alltoallv(out, numg, startg, MPI_COMPLEX, in, numr, startr, MPI_COMPLEX, this->pool_world);
    ModuleBase::GlobalFunc::ZEROS(out, this->nrxx);
    //change (nplane,nstot) to (nplane fftnxy)
	for (int istot = 0;istot < nstot; ++istot)
	{
		int ixy = this->istot2ixy[istot];
        //int ixy = (ixy / fftny)*ny + ixy % fftny;
        outp = &out[ixy * nplane];
        inp = &in[istot * nplane];
		for (int iz = 0; iz < nplane; ++iz)
		{
			outp[iz] = inp[iz];
		}
	}

#endif
    ModuleBase::timer::tick(this->classname, "gathers_scatterp");
    return;
}



}