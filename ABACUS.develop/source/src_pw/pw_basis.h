//==========================================================
// AUTHOR : Lixin He, Fang Wei, Mohan Chen
// DATA : 2006-11 ~ 2008-11
//==========================================================
#ifndef PLANEWAVE_H
#define PLANEWAVE_H

#include "../src_global/complexmatrix.h"
#include "../src_global/vector3.h"
#include "unitcell.h"
#include "klist.h"
#include "../src_parallel/ft.h"
using namespace std;

class PW_Basis
{
	public:

    PW_Basis();
    ~PW_Basis();

	const UnitCell *Ucell; // pointer for UnitCell
	const kvect *Klist; // pointer for K-point list
	FFT FFT_chg; // FFT operations for charge density
	FFT FFT_wfc; // FFT operations for wave functions

    void set
    (
        const bool &gamma_only_in,
        const double &ecutwfc_in,
        const double &ecutrho_in,
        const int &nx_in,
        const int &ny_in,
        const int &nz_in,
        const int &ncx_in,
        const int &ncy_in,
        const int &ncz_in,
		const int &bx_in,
		const int &by_in,
		const int &bz_in
    );

    //===============================================
    // Part 1: Kinetic energy cutoff for plane waves
    //===============================================
    bool gamma_only;	// only half G are used.
    double ecutwfc;		// Energy cutoff for wavefunctions.
    double ecutrho; 	// Energy cutoff for charge/potential.

    // ggpsi=2*Ecut*(lat0*lat0)/(4.0*PI*PI);
    // in NCPP,  ggchg=ggfft; in USPP (double grids), ggchg=4*ggfft
    double ggpsi; 		// planewave cut off for the wavefunctions, unit (2*PI/lat0)^2
    double ggwfc;		// cutoff for wave function, ggwfc=4*ggpsi
    double ggwfc2;		// ggwav=wfact*ggpsi, default value: wfact=4.0
    double wfac;		// weighting factor from ggpsi to ggwfc2
    double ggchg;   	// G^2 cutoff for supporting charge density,

    //===============================================
    // Part 2: FFT dimensions in real space
    //===============================================

	//-----------------------------------------------
	// FFT dimensions for wave functions.
	int nx, ny, nz, nxyz;
	// FFT dimensions for charge/potential.
	int ncx, ncy, ncz, ncxyz;
	// nczp: number of x-y planes in each processor
	// nczp_start: starting nczp index in each processor
	int nczp, nczp_start;
	// nrxx: for parallel,local FFT grid dimension
    int nrxx;
	// nrxx_start: starting nrxx in each processor 
    int nrxx_start;

	//-----------------------------------------------
	// used in LCAO algorithm: grid integrals 
	// combine [bx,by,bz] FFT grids into a big one
	// typical values are bx=2, by=2, bz=2
	// bxyz = bx*by*bz, typical value is 8
    int bx, by, bz, bxyz;
	// nbx=ncx/bx, nby=ncy/by, nbz=ncz/bz, 
	// nbxyz=nbx*nby*nbz
	int nbx, nby, nbz, nbxyz;
	// nbzp: number of x-y planes in each processor 
	// nbzp_start: starting nbzp index in each processor
	int nbzp, nbzp_start;
	// nbxx is local dimension of nbxyz
	// nbxx=nbx*nby*nbzp
	int nbxx;


    //===============================================
    // Part 3: G vectors in reciprocal FFT box
    //===============================================

	// structure factor (ntype, ngmax)
    ComplexMatrix strucFac;

	// (= ngmax) / num. of G vectors within ggfft
    int ngmw;

	// mohan add 2008-3-25
	// (=*ind_FFT) /1D G vector -> for wave function FFT
    int *ig2fftw;
    //ig2fftw= new int [ngmw] for wave functions
    //ig2fftw[ng] -> the coordinates in FFT box

    //===============================================
    // Part 4: Charge & Potential (G vectors,FFT box)
    //===============================================

    // G vectors, FFT box related to charge density and potentials etc.
    // FFT grid for charge density, potentials etc.
    // this grid equal to the wave function grid


	// new, num. of G vectors within ggchg
    int ngmc;

	// mohan add 2008-4-2
    int ngmc_g;

	// mohan add 2008-4-2
    int *ig2fftc;
    //ig2fftc= new int [ngmc] for wave functions
    //ig2fftc[ng] -> the coordinates in FFT box

    Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc],
    Vector3<double> *gdirect_global;	//(= *G1d) ; // ig = new Vector igc[ngmc],
    // store the 3D G vector coordinates for charge density grid;
    //ig.x, ig.y, ig.z should be integers !!!
    // G vectors are in order of increasing G^2
    // can be shared by charge density/potential and wave functions.

    Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    Vector3<double> *gcar_global;   	//G vectors in cartesian corrdinate
    //g=ig*G ?? HLX (05-26-06): need to check if this is ok!
    //tau is also defined in Cartesian coordintes unit lat0

    double *gg;         			// modulus (G^2) of G vectors
    double *gg_global;    			// modulus (G^2) of G vectors
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)

    //==================
    // void set_nggm()
    //==================
    int nggm;          				// =ngl(another dft code);

    // number of |G| shells (i.e., how many different |G|s)
    double *ggs;         			// store the |G|^2 for each shell
    int *ig2ngg;        			//=*igtong1 (another dft code)
    // G -> correspondence shells of G
    // end charge grids

    int gstart;        				// first nonzero g vector

    ComplexMatrix eigts1;   		//
    ComplexMatrix eigts2;   		//the phases e^{-iG*tau_s}
    ComplexMatrix eigts3;   		//

    int *ig1;        				//
    int *ig2;        				// the indices of G components
    int *ig3;        				//

    double *gg_global0; //LiuXh add 20180515
    int *cutgg_num_table; //LiuXh add 20180515
    int ggchg_time_global; //LiuXh add 20180515


    void gen_pw(ofstream &log, const UnitCell &Ucell_in, const kvect &Klist_in);

    void setup_structure_factor(void); 		// Calculate structur factors

    void update_gvectors(ofstream &log, const UnitCell &Ucell_in); //LiuXh add 20180515

	private:

    void setup_gg(void);

    void setup_FFT_dimension(void); 		// set up FFT dimensions

#ifdef __MPI
    void divide_fft_grid(void);

    void get_MPI_GVectors(void);

    void columns_and_pw_distribution_2(void);
#else

    void get_GVectors(void);

#endif

    void get_nggm(const int ngmc_local);

};
#endif //PlaneWave class
