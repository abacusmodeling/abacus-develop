#ifndef PWBASIS_H
#define PWBASIS_H

#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/vector3.h"
#include <assert.h>
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
        ModuleBase::Matrix3 G_in, // reciprocal lattice vector (2pi*inv(R) )
        double gridecut
    );
    //Init the grids for FFT
    void initgrids(
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        ModuleBase::Matrix3 G_in, // reciprocal lattice vector (2pi*inv(R) )
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

    //distribute plane waves to different processors
    void distribute_g();

    //distribute real-space grids to different processors
    void distribute_r();
//===============================================
// distribution maps
//===============================================
public:
    //reciprocal-space
    int *ig2isz; // map ig to (is, iz).
    int *istot2ixy; // istot2ixy[is]: ix + iy * nx of is^th stick among all sticks.
    int *is2ixy; // is2ixy[is]: ix + iy * nx of is^th stick among sticks on current proc.
    int *ixy2ip; // store the ip of proc which contains stick on (x, y).
    int *startis; // startis[ip]: starting is stick in the ip^th proc.
    int *nst_per; // number of sticks on each core.
    int nst; //num. of sticks in current proc.
    int npw; //num. of plane waves in current proc.
    //real space
    int nrxx;
    int *startz; //startz[ip]: starting z plane in the ip-th proc. in current POOL_WORLD 
	int *numz; //numz[ip]: num. of z planes in the ip-th proc. in current POOL_WORLD


   

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
    void divide_pw(
        const int tot_npw,                          // total number of planewaves.
        double* gg_global,                          // the modulus of all planewaves.
        ModuleBase::Vector3<double>*gdirect_global, // direct coordinates of planewaves.
        double* gg2D,                               // the i^th row contains the modulus of planewaves that belong to the i^th core.
        ModuleBase::Vector3<double>*gdirect2D       // the i^th row contains the direct coordinates of planewaves that belong to the i^th core.
    );
    void get_ig2fft_is2ir(    
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j,          // y or y + ny (if y < 0) of stick.
        int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
        int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
        int* is2ip,         // ip of core containing is^th stick, map is to ip.
        int* nst_per       // number of sticks on each core.
    );


//===============================================
// Part 2: FFT dimensions in real space
//===============================================
public:
	// FFT dimensions for wave functions.
	int nx, ny, nz, nxyz;
	// FFT dimensions for charge/potential.
    
    int seed;
    
    
    ModuleBase::Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc]
    ModuleBase::Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    double *gg;       	// modulus (G^2) of G vectors [ngmc]
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)
	// gg_global dimension: [cutgg_num_now] (save memory skill is used)


};
#endif //PlaneWave class
