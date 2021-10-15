#ifndef PWBASIS_H
#define PWBASIS_H

#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/vector3.h"
#include "../src_pw/klist.h"
#include <assert.h>
//
//A class which can convert a function of "r" to the corresponding linear 
// superposition of plane waves (real space to reciprocal space)
// or convert a linear superposition of plane waves to the function 
// of "r" (reciprocal to real).
//plane waves: <r|g,k>=1/sqrt(V)*exp(i(k+g)r)
//
class PW_Basis
{

public:
    PW_Basis();
    ~PW_Basis();

    //Init the grids for FFT
    void initgrids(
        ModuleBase::Matrix3 latvec_in; // Unitcell lattice vectors
        ModuleBase::Matrix3 G_in; // reciprocal lattice vector (2pi*inv(R) )
        double gridecut
    );
    //Init the grids for FFT
    void initgrids(
        ModuleBase::Matrix3 latvec_in; // Unitcell lattice vectors
        ModuleBase::Matrix3 G_in; // reciprocal lattice vector (2pi*inv(R) )
        int nx_in, int ny_in, int nz_in
    );

    //Init some parameters
    void initparameters(
        bool gamma_only_in,
        double ecut_in,
        int nk_in, //number of k points in this pool
        ModuleBase::Vector3<double> *kvec_d, // Direct coordinates of k points
        int poolnproc_in, // Number of processors in this pool
        int poolrank_in, // Rank in this pool
        int distribution_type_in,
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
    int *ig2fft; // dimension: [ngmw]
    int *is2ir;
    int nst; //num. of sticks in current proc.
    int npw; //num. of plane waves in current proc.
    //real space
    int nrxx;
    int *startz; //startz[ip]: starting z plane in the ip-th proc. in current POOL_WORLD 
	int *numz; //numz[ip]: num. of z planes in the ip-th proc. in current POOL_WORLD


   

private:
    bool gamma_only;	// only half g are used.
    double ggecut;    //Energy cut off for g^2 
    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors, unit in bohr
    ModuleBase::Matrix3 G; // reciprocal lattice vector, unit in 
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    int distribution_type;
    int poolnproc;
    int poolrank;
    int nks;
    ModuleBase::Vector3<double> *kvec_d;

    //distribute plane waves to different processors
    void distribution_method1();


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
