#ifndef PWBASIS_H
#define PWBASIS_H

#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/vector3.h"
#include <complex>
#include "fft.h"

namespace ModulePW
{

/**
 * @brief A class which can convert a function of "r" to the corresponding linear 
 * superposition of plane waves (real space to reciprocal space)
 * or convert a linear superposition of plane waves to the function 
 * of "r" (reciprocal to real).
 * @author qianrui, Sunliang on 2021-10-15
 * @details
 * Math:
 * plane waves: <r|g>=1/sqrt(V) * exp(igr)
 * f(r) = 1/sqrt(V) * \sum_g{c(g)*exp(igr)}
 * c(g) = \int f(r)*exp(-igr) dr
 * USAGE：
 * ModulePW::PW_Basis pwtest;
 * 1. setup FFT grids for PW_Basis
 * pwtest.initgrids(lat0,latvec,gridecut,nproc_in_pool,rank_in_pool);
 * pwtest.initgrids(lat0,latvec,N1,N2,N3,nproc_in_pool,rank_in_pool); 
 * //double lat0：unit length, (unit: bohr)
 * //ModuleBase::Matrix3 latvec：lattice vector, (unit: lat0), e.g. ModuleBase::Matrix3 latvec(1, 1, 0, 0, 2, 0, 0, 0, 2);
 * //double gridecut：cutoff energy to generate FFT grids, (unit: Ry)
 * //int N1,N2,N3: FFT grids
 * 2. init parameters
 * pwtest.initparameters(gamma_only,ggecut,dividemthd);
 * //bool gamma_only: if use gamma_only
 * //double ggecut: cutoff kinetic energy for planewaves,(unit in Ry) G^2 < ggecut
 * //int dividemthd: method to divide planewaves to different cores
 * 3. Setup transforms from real space to reciprocal space or from reciprocal space to real space.
 * pwtest.setuptransform(); 
 * pwtest.recip2real(rhog,rhor); //rhog to rhor
 * pwtest.real2recip(rhor,rhog); //rhor to rhog
 * 4. Generate the wave vector for planewaves
 * pwtest.collect_local_pw(); 
 * //then we can use pwtest.gg, pwtest.gdirect, pwtest.gcar, (unit in lat0^-1 or lat0^-2)
 * 
 */
class PW_Basis
{

public:
    PW_Basis();
    ~PW_Basis();

    //Init the grids for FFT
    void initgrids(
        const double lat0_in, //unit length (unit in bohr)
        const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors (unit in lat0) 
        const double gridecut, //unit in Ry, ecut to set up grids
        const int poolnproc_in, // Number of processors in this pool
        const int poolrank_in // Rank in this pool
    );
    //Init the grids for FFT
    void initgrids(
        const double lat0_in,
        const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        const int nx_in, int bigny_in, int nz_in,
        const int poolnproc_in, // Number of processors in this pool
        const int poolrank_in // Rank in this pool
    );

    //Init some parameters
    void initparameters(
        const bool gamma_only_in,
        const double pwecut_in, //unit in Ry, ecut to decides plane waves
        const int distribution_type_in = 1
    );

//===============================================
//                 distribution maps
//===============================================
public:
    //reciprocal-space
    // only on first proc.
    int *startnsz_per;//useless // startnsz_per[ip]: starting is * nz stick in the ip^th proc.
    int *nstnz_per;//useless // nz * nst(number of sticks) on each core.
    int *nst_per;// nst on each core
    // on all proc.
    int *ig2isz; // map ig to (is, iz).
    int *istot2bigixy; // istot2bigixy[is]: iy + ix * bigny of is^th stick among all sticks.
    int *ixy2istot; //useless // ixy2istot[iy + ix * ny]: is of stick on (ix, iy) among all sticks.
    int *is2ixy; // is2ixy[is]: iy + ix * bigny of is^th stick among sticks on current proc.
    int *ixy2ip; // useless// ixy2ip[iy + ix * ny]: ip of proc which contains stick on (ix, iy).
    int nst; //num. of sticks in current proc.
    int nstnz; // nst * nz
    int nstot; //num. of sticks in total.
    int npw; //num. of plane waves in current proc.
    //real space
    int nrxx; //num. of real space grids
    int *startz; //startz[ip]: starting z plane in the ip-th proc. in current POOL_WORLD 
	int *numz; //numz[ip]: num. of z planes in the ip-th proc. in current POOL_WORLD
    int *numg; //numg[ip] :  nst_per[poolrank] * numz[ip] 
    int *numr; //numr[ip] :  numz[poolrank] * nst_per[ip]
    int *startg;  // startg[ip] = numg[ip-1] + startg[ip-1]
    int *startr;  // startr[ip] = numr[ip-1] + startr[ip-1]
    int nplane; //num. of planes in current proc.

    ModuleBase::Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc]
    ModuleBase::Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    double *gg;       	// modulus (G^2) of G vectors [ngmc]
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)
	// gg_global dimension: [cutgg_num_now] (save memory skill is used)

    //distribute plane waves and grids and set up fft
    void setuptransform();

protected:
    //distribute plane waves to different processors
    void distribute_g();

    //distribute real-space grids to different processors
    void distribute_r();

    //prepare for MPI_Alltoall
    void getstartgr();


public:
    //prepare for transforms between real and reciprocal spaces
    void collect_local_pw();

    // void collect_tot_pw(
    //     double* gg_global,
    //     ModuleBase::Vector3<double> *gdirect_global,
    //     ModuleBase::Vector3<double> *gcar_global
    // ); 
   

public:
    bool gamma_only;	// only half g are used.
    double ggecut;    //Energy cut off for g^2/2, unit in 1/lat0^2, ggecut=ecutwfc(Ry)*lat0^2/4pi^2
    double lat0;     //unit length for lattice, unit in bohr
    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors, unit in lat0
    ModuleBase::Matrix3 G; // reciprocal lattice vector, unit in 1/lat0
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    int distribution_type;
    int poolnproc;
    int poolrank;
    //distribute plane waves to different processors

protected:
    //method 1: first consider number of plane waves
    void distribution_method1(); 
    // Distribute sticks to cores in method 1.
    void divide_sticks_1(
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j,          // y or y + ny (if y < 0) of stick.
        int* st_length,     // the stick on (ix, iy) consists of st_length[ix*ny+iy] planewaves.
        int* npw_per       // number of planewaves on each core.
    );

    //method 2: first consider number of sticks
    void distribution_method2();
    // Distribute sticks to cores in method 2.
    void divide_sticks_2();
   
    //Count the total number of planewaves (tot_npw) and sticks (this->nstot) (in distributeg method1 and method2)
    void count_pw_st(
        int &tot_npw,     // total number of planewaves.
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
    );

    //get ig2isz and is2ixy
    void get_ig2isz_is2ixy(
        int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
        int* st_length     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    );

    //Collect the x, y indexs, length of the sticks (in distributeg method1)
    void collect_st(
        int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
        int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
        int* st_i,                                      // x or x + nx (if x < 0) of stick.
        int* st_j,                                      // y or y + ny (if y < 0) of stick.
        int* st_length                                 // number of planewaves in stick, stored in 1d array with tot_nst elements.
    );

    //get istot2bigixy
    void get_istot2bigixy(
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j          // y or y + ny (if y < 0) of stick.
    );

    //Create the maps from ixy to (in method 2)
    void create_maps(
        int* st_length2D,  // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
        int* npw_per       // number of planewaves on each core.
    );

//===============================================
//                  FFT
//===============================================
public:
	// FFT dimensions for wave functions.
	int nx, ny, nz, nxyz, nxy;
    int bigny, bignxyz, bignxy; // Gamma_only: ny = int(bigny/2)-1 , others: ny = bigny
    int liy,riy;// liy: the left edge of the pw ball; riy: the right edge of the pw ball
    int maxgrids; // max between nz * ns and bignxy * nplane
    FFT ft;

    void real2recip(double * in, std::complex<double> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(std::complex<double> * in, std::complex<double> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(std::complex<double> * in, double *out); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(std::complex<double> * in, std::complex<double> * out); //in:(nz, ns)  ; out(nplane,nx*ny)

#ifdef __MIX_PRECISION
    void real2recip(float * in, std::complex<float> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(std::complex<float> * in, std::complex<float> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(std::complex<float> * in, float *out); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(std::complex<float> * in, std::complex<float> * out); //in:(nz, ns)  ; out(nplane,nx*ny)
#endif
protected:
    //gather planes and scatter sticks of all processors
    template<typename T>
    void gatherp_scatters(std::complex<T> *in, std::complex<T> *out); 

    //gather sticks of and scatter planes of all processors
    template<typename T>
    void gathers_scatterp(std::complex<T> *in, std::complex<T> *out); 
};

}
#endif //PlaneWave 
