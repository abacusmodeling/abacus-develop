#ifndef PWBASIS_H
#define PWBASIS_H

#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/vector3.h"
#include <complex>
#include "fft.h"

//
//A class which can convert a function of "r" to the corresponding linear 
// superposition of plane waves (real space to reciprocal space)
// or convert a linear superposition of plane waves to the function 
// of "r" (reciprocal to real).
//plane waves: <r|g>=1/sqrt(V)*exp(igr)
//
class PW_Basis
{

public:
    PW_Basis();
    ~PW_Basis();

    //Init the grids for FFT
    void initgrids(
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        double gridecut
    );
    //Init the grids for FFT
    void initgrids(
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        int nx_in, int ny_in, int nz_in
    );

    //Init some parameters
    void initparameters(
        bool gamma_only_in,
        double ggecut_in,
        int poolnproc_in, // Number of processors in this pool
        int poolrank_in, // Rank in this pool
        int distribution_type_in
    );

//===============================================
//                 distribution maps
//===============================================
public:
    //reciprocal-space
    // only on first proc.
    int *startnsz_per; // startnsz_per[ip]: starting is * nz stick in the ip^th proc.
    int *nstnz_per; // nz * nst(number of sticks) on each core.
    // on all proc.
    int *ig2isz; // map ig to (is, iz).
    int *istot2ixy; // istot2ixy[is]: ix + iy * nx of is^th stick among all sticks.
    int *ixy2istot; // ixy2istot[ix + iy * nx]: is of stick on (ix, iy) among all sticks.
    int *is2ixy; // is2ixy[is]: ix + iy * nx of is^th stick among sticks on current proc.
    int *ixy2ip; // ixy2ip[ix + iy * nx]: ip of proc which contains stick on (ix, iy).
    int nst; //num. of sticks in current proc.
    int nstot; //num. of sticks in total.
    int npw; //num. of plane waves in current proc.
    //real space
    int nrxx; //num. of real space grids
    int *startz; //startz[ip]: starting z plane in the ip-th proc. in current POOL_WORLD 
	int *numz; //numz[ip]: num. of z planes in the ip-th proc. in current POOL_WORLD

    ModuleBase::Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc]
    ModuleBase::Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    double *gg;       	// modulus (G^2) of G vectors [ngmc]
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)
	// gg_global dimension: [cutgg_num_now] (save memory skill is used)

    //distribute plane waves to different processors
    void distribute_g();

    //distribute real-space grids to different processors
    void distribute_r();

    //distribute plane waves to different processors
    void distribution_method1();


   

private:
    bool gamma_only;	// only half g are used.
    double ggecut;    //Energy cut off for g^2/2 
    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors, unit in bohr
    ModuleBase::Matrix3 G; // reciprocal lattice vector, unit in 
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    int distribution_type;
    int poolnproc;
    int poolrank;
   

    //distribute plane waves to different processors
    void distribution_method1();
    void count_pw_st(
        int &tot_npw,     // total number of planewaves.
        int &tot_nst,     // total number of sticks.
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
    );
    void collect_pw_st(
        const int tot_npw,                              // total number of planewaves.
        const int tot_nst,                              // total number of sticks.
        int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
        int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
        double* gg_global,                              // the modulus of all planewaves.
        ModuleBase::Vector3<double> *gdirect_global,    // direct coordinates of planewaves.
        int* st_i,                                      // x or x + nx (if x < 0) of stick.
        int* st_j,                                      // y or y + ny (if y < 0) of stick.
        int* st_length,                                 // number of planewaves in stick, stored in 1d array with tot_nst elements.
        int* st_bottom                                  // minimum z of stick, stored in 1d array with tot_nst elements.
    );
    void divide_sticks(
        const int tot_npw,  // total number of planewaves.
        const int tot_nst,  // total number of sticks.
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j,          // y or y + ny (if y < 0) of stick.
        int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
        int* npw_per,       // number of planewaves on each core.
        int* nst_per,       // number of sticks on each core.
        int* is2ip         // ip of core containing is^th stick, map is to ip.         
    );
    void get_istot2ixy(
        const int tot_nst,  // total number of sticks.
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j,          // y or y + ny (if y < 0) of stick.
        int* is2ip          // ip of core containing is^th stick, map is to ip.
    );
    void divide_pw(
        const int tot_npw,                          // total number of planewaves.
        double* gg_global,                          // the modulus of all planewaves.
        ModuleBase::Vector3<double>*gdirect_global, // direct coordinates of planewaves.
        double** gg2D,                               // the i^th row contains the modulus of planewaves that belong to the i^th core.
        ModuleBase::Vector3<double>**gdirect2D       // the i^th row contains the direct coordinates of planewaves that belong to the i^th core.
    );
    void get_ig2isz_is2ixy(
        const int tot_nst,  // total number of sticks.
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j,          // y or y + ny (if y < 0) of stick.
        int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
        int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
        int* is2ip          // ip of core containing is^th stick, map is to ip.
    );

//===============================================
//                  FFT
//===============================================
public:
	// FFT dimensions for wave functions.
	int nx, ny, nz, nxyz, nxy;
    FFT ft;

    void real2recip(double * in, complex<double> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(complex<double> * in, complex<double> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(complex<double> * in, double *out); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(complex<double> * in, complex<double> * out); //in:(nz, ns)  ; out(nplane,nx*ny)

    void gatherp_scatters(complex<double> *in, complex<double> *out); //gather planes and scatter sticks of all processors
    void gathers_scatterp(complex<double> *in, complex<double> *out); //gather sticks of and scatter planes of all processors
    void gatherp_scatters_gamma(complex<double> *in, complex<double> *out); //gather planes and scatter sticks of all processors, used when gamma_only
    void gathers_scatterp_gamma(complex<double> *in, complex<double> *out); //gather sticks of and scatter planes of all processors, used when gamma only

    
    
   


};
#endif //PlaneWave class
