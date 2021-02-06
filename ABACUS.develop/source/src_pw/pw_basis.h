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

	// pointer for UnitCell
	const UnitCell *Ucell;

	// pointer for K-point list
	const kvect *Klist;

	// FFT grid for charge density
	FFT FFT_chg;

	// FFT grid for wave functions
	FFT FFT_wfc;

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

    //======================
    // Part1:Energy cutoff
    //======================
    bool gamma_only;				// only half G are used.
    double wfac;					// weighting factor
    double ecutwfc;					// Energy cutoff for wavefunctions.
    double ecutrho; 				// Energy cutoff for charge/potential.

    // ggpsi=2*Ecut*(lat0*lat0)/(4.0*PI*PI);
    // in NCPP,  ggchg=ggfft
    // in USPP (double grids), ggchg=4*ggfft
    double ggpsi; 					// planewave cut off for the wavefunctions, unit (2*PI/lat0)^2
    double ggwfc;					// ggwav;(=G2max);/ G^2 cutoff for wave function FFT box, unit (2*PI/a0)^2,
    double ggwfc2;					// ggwav=wfact*ggpsi, default value: wfact=4.0
    double ggchg;   				// G^2 cutoff for supporting charge density,

    //===============
    // FFT_dimension
    //===============
    int bx, by, bz, bxyz;			// subset of each grid of FFT box, for example 2 2 2
	int nbx, nby, nbz, nbxyz;		//
	int nbzp, nbzp_start;			// start ncz in each process, mohan add 2009-11-09
	int nbxx, nbxx_start;

	int nx, ny, nz, nxyz;			// fft dimensions for wave functions.
	int ncx, ncy, ncz, ncxyz;		// fft dimensions for charge/potential.
	int nczp, nczp_start;			// start ncz in each process, mohan add 2009-11-09

    ComplexMatrix strucFac;			// StrucFac (ntype,ngmax)

    //=========================================
    // Part2: Wave function(G vectors,FFT box)
    //=========================================
    int ngmw;						//(= ngmax) / num. of G vectors within ggfft
    /*** mohan add 2008-3-25 ***/
    int *ig2fftw;					//(=*ind_FFT) /1D G vector -> for wave function FFT
    //ig2fftw= new int [ngmw] for wave functions
    //ig2fftw[ng] -> the coordinates in FFT box

    //===============================================
    // Part3: Charge & Potential (G vectors,FFT box)
    //===============================================
    // G vectors, FFT box related to charge density and potentials etc.
    // fft grid for charge density, potentials etc.
    // this grid equal to the wave function grid
    int nrxx;						//for parallel,local FFT grid dimension
    int nrxx_start;                 //for parallel,local FFT grid dimension
    int ngmc;       				// new, num. of G vectors within ggchg,
    int ngmc_g;						// mohan add 2008-4-2
    int *ig2fftc;					//(=*ind_FFT) //1D G vector -> charge/potential FFT
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

    bool doublegrid;				// .TRUE. if we use a double grid; (for USPP)
    // .FALSE., if NCPP.

    void gen_pw(ofstream &log, const UnitCell &Ucell_in, const kvect &Klist_in);

    void setup_structure_factor(); 		// Calculate structur factors

    void update_gvectors(ofstream &log, const UnitCell &Ucell_in); //LiuXh add 20180515

	private:

    void setup_gg();

    void setup_FFT_dimension(); 		// set up FFT dimensions

#ifdef __MPI
    void divide_fft_grid();

    void get_MPI_GVectors();

    void columns_and_pw_distribution_2();
#else

    void get_GVectors();

#endif

    void get_nggm(const int ngmc_local);

};
#endif //PlaneWave class
