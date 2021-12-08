#include "pw_basis.h"
#include "../module_base/global_function.h"
#ifdef __MPI
#include "mpi.h"
#include "../src_parallel/parallel_global.h"
#endif
namespace ModulePW
{

//in: (nplane,nxy) out:(nz,nst)
void PW_Basis:: gatherp_scatters(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = this->is2ixy[is];
            int bigixy = (ixy / ny)*bigny + ixy % ny;
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                out[is*nz+iz] = in[bigixy*nz+iz];
            }
        }
        return;
    }
#ifdef __MPI
    std::complex<double> * tmp = NULL;
    if(this->poolrank == 0) tmp = new std::complex<double> [this->nz * this->nstot];
    
    //gather planes of different processors
    for(int ixy = 0 ; ixy < this->nxy ; ++ixy)
    {
        if(this->ixy2ip[ixy] == -1) continue;
        int istot = 0;
        if(this->poolrank == 0) istot = this->ixy2istot[ixy];
        int bigixy = (ixy / ny)*bigny + ixy % ny;
        MPI_Gatherv(&in[bigixy*this->nplane], this->nplane, mpicomplex, &tmp[istot*this->nz], 
                    this->numz,this->startz,mpicomplex,0,POOL_WORLD);
    }
    
    //scatter sticks to different processors
    MPI_Scatterv(tmp, this->nstnz_per, this->startnsz_per,mpicomplex,out,
                    this->nstnz,mpicomplex,0, POOL_WORLD); 
    
    if(tmp!=NULL) delete[] tmp;
#endif
    return;
}


//in: (nz,nst) out:(nplane,nxy)
void PW_Basis:: gathers_scatterp(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) //In this case nrxx=nx*ny*nz, nst = nstot, 
    {
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] = 0.0;
        }
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = is2ixy[is];
            int bigixy = (ixy / ny)*bigny + ixy % ny;
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                out[bigixy*nz+iz] = in[is*nz+iz];
            }
        }
        return;
    }
#ifdef __MPI
    if(this->poolnproc == 1) return;
    std::complex<double> * tmp = NULL;
    if(this->poolrank == 0) tmp = new std::complex<double> [this->nz * this->nstot];

    //scatter sticks to different processors
    MPI_Gatherv(in, this->nstnz, mpicomplex, tmp,
                    this->nstnz_per, this->startnsz_per, mpicomplex, 0, POOL_WORLD);

    for(int ir = 0 ; ir < this->nrxx ; ++ir) out[ir] = 0.0;
    for(int istot = 0 ; istot < this->nstot ; ++istot)
    {
        int ixy = this->istot2ixy[istot];
        int bigixy = (ixy / ny)*bigny + ixy % ny;
        MPI_Scatterv(&tmp[istot*this->nz], this->numz,this->startz, mpicomplex, &out[bigixy*this->nplane], 
                    this->nplane,mpicomplex,0,POOL_WORLD);
    }
    if(tmp!=NULL) delete[] tmp;
#endif
    return;
}

//in: (nplane,nxy) out:(nz,nst)
void PW_Basis:: gatherp_scatters2(complex<double> *in, complex<double> *out)
{
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = this->is2ixy[is];
            int bigixy = (ixy / ny)*bigny + ixy % ny;
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                out[is*nz+iz] = in[bigixy*nz+iz];
            }
        }
        return;
    }
#ifdef __MPI
    std::complex<double>*tmp = new std::complex<double> [this->maxgrids];
    //change (nplane nxy) to (nplane,nstot)
    // Hence, we can send them at one time.  
	for (int istot = 0;istot < nstot; ++istot)
	{
		int ixy = this->istot2ixy[istot];
        int bigixy = (ixy / ny)*bigny + ixy % ny;
		for (int iz = 0; iz < nplane; ++iz)
		{
			tmp[istot*nplane+iz] = in[bigixy*nplane+iz];
		}
	}

    //exchange data
    //(nplane,nstot) to (numz[ip],ns, poolnproc)
	MPI_Alltoallv(tmp, numr, startr, mpicomplex, out, numg, startg, mpicomplex, POOL_WORLD);

    // change (nz,ns) to (numz[ip],ns, poolnproc)
    for (int ip = 0; ip < this->poolnproc ;++ip)
	{
        int nzip = this->numz[ip];
		for (int is = 0; is < this->nst; ++is)
		{
			for (int izip = 0; izip < nzip; ++izip)
			{
				tmp[ is * nz + startz[ip] + izip] = out[startg[ip] + is*nzip + izip];
			}
		}
	}

    for(int isiz = 0 ; isiz < this->nz * this->nst ; ++isiz)
    {
        out[isiz] = tmp[isiz];
    }
    
   
#endif
    return;
}

//in: (nz,nst) out:(nplane,nxy)
//in[] will be changed
void PW_Basis:: gathers_scatterp2(complex<double> *in, complex<double> *out)
{
    
    if(this->poolnproc == 1) //In this case nrxx=nx*ny*nz, nst = nstot, 
    {
        ModuleBase::GlobalFunc::ZEROS(out, this->nrxx);
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = is2ixy[is];
            int bigixy = (ixy / ny)*bigny + ixy % ny;
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                out[bigixy*nz+iz] = in[is*nz+iz];
            }
        }
        return;
    }
#ifdef __MPI
    // std::complex<double>*tmp = new std::complex<double> [this->maxgrids];
    // for(int isiz = 0 ; isiz < this->nz * this->nst ; ++isiz)
    // {
    //     tmp[isiz] = in[isiz];
    // }

    // change (nz,ns) to (numz[ip],ns, poolnproc)
    // Hence, we can send them at one time.  
    for (int ip = 0; ip < this->poolnproc ;++ip)
	{
        int nzip = this->numz[ip];
		for (int is = 0; is < this->nst; ++is)
		{
			for (int izip = 0; izip < nzip; ++izip)
			{
				out[startg[ip] + is*nzip + izip] = in[ is * nz + startz[ip] + izip];
			}
		}
	}

	//exchange data
    //(numz[ip],ns, poolnproc) to (nplane,nstot)
	MPI_Alltoallv(out, numg, startg, mpicomplex, in, numr, startr, mpicomplex, POOL_WORLD);

    ModuleBase::GlobalFunc::ZEROS(out, this->nrxx);
    //change (nplane,nstot) to (nplane nxy)
	for (int istot = 0;istot < nstot; ++istot)
	{
		int ixy = this->istot2ixy[istot];
        int bigixy = (ixy / ny)*bigny + ixy % ny;
		for (int iz = 0; iz < nplane; ++iz)
		{
			out[bigixy*nplane+iz] = in[istot*nplane+iz];
		}
	}

    // delete[] tmp;
#endif
    return;
}



}