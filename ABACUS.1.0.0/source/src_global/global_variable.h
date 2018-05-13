//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef GLOBAL_VARIABLE_H
#define GLOBAL_VARIABLE_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
//==========================================================
// EXPLAIN : Basic Global Variables
//==========================================================

extern int 		NBANDS;				// 1
extern int		NBANDS_ISTATE;		// 1.05 // mohan add 2011-03-22
extern int 		NLOCAL;				// 1.1 // mohan add 2009-05-29

extern string	CALCULATION;		// 2 "scf";"nscf" ;"symmetry"
extern bool 	BERRY_PHASE;		// 3 berry phase calculation
extern bool 	LDA_PLUS_U;			// 4 LDA + U method
extern int		EFIELD;				// 5 add electric field
extern int		BFIELD;				// 6 add magnetic field
extern int 		DIPOLE;				// 7 add dipole correction


extern string   DFT_FUNCTIONAL;		// 6.5 change the DFT functional from input file.
extern int 		NSPIN;				// 7
extern bool		TWO_EFERMI; 		// 7.5 mohan add 2011-04-03, two fermi energy, exist if magnetization is fixed.
extern int 		CURRENT_SPIN;		// 8

extern int 		VNA;				// 8.01 // mohan add 2011-05-27
extern int		GRID_SPEED;			// 8.02 // mohan add 2012-03-29
extern int 		FORCE;				// 8.1
extern double	FORCE_THR;			// 8.2
extern bool 	STRESS;				// 8.25 calcualte the stress
extern string	MOVE_IONS;			// 8.26
extern string	OUT_LEVEL;			// 8.27 

extern int		NSTEP;				// 8.3
extern int 		NITER;				// 8.4

extern bool 	SYMMETRY;		// 9
extern bool		MLWF_FLAG;		// 9.1 mohan add 2010-01-26
//extern int 		LOCAL_BASIS; xiaohui modify 2013-09-01		// 10
//extern int	 	LINEAR_SCALING;	xiaohui modify 2013-09-01	// 11 // -1: test 0: pw 1:lcao 2:O(N)
extern string		BASIS_TYPE; //xiaohui add 2013-09-01
extern string		KS_SOLVER; //xiaohui add 2013-09-01
extern double	SEARCH_RADIUS;		// 11.1 // mohan add 2011-03-10
extern bool		SEARCH_PBC;			// 11.2 // mohan add 2011-03-10
extern bool		SPARSE_MATRIX;		// 11.3 // mohan add 2009-03-13
extern int 		ATOM_DISTRIBUTION;	// 11.4 // mohan add 2010-06-28

//added by zhengdy-soc
extern bool     NONCOLIN;
extern bool     LSPINORB;
extern bool     DOMAG;
extern int      NPOL;

// diagonalization (5)
//extern string 	DIAGO_TYPE; xiaohui modify 2013-09-01			// 12 "cg","davidson","fs","hpseps"
extern int		DIAGO_PROC;			// 12.1 number of processors used to diag.
extern int 		DIAGO_CG_MAXITER;	// 13
extern int		DIAGO_CG_PREC;		// 13.1
extern int 		DIAGO_DAVID_NDIM;	// 14
extern double 	ETHR;				// 15 ethr
extern double	FS_REF_ENERGY;		// 16
extern int		NB2D;				// 16.5 dividsion of 2D_matrix.

extern double 	DRHO2;				// 17

extern string	RESTART_MODE;		// 18

extern double DQ;	// 19 mohan add 2009-09-10
extern int NQX;	// 20 mohan add 2009-09-10

extern int NURSE; // 21 mohan add 2010-09-10
extern bool COLOUR; // mohan add 2011-04-26
extern bool GAMMA_ONLY_LOCAL; // 22 : mohan add 2010-10-20
extern bool GAMMA_ONLY_PW; // mohan add 2012-06-05

extern int T_IN_H; // 23, calculate T in H or not.
extern int VL_IN_H; // 24, calculate Vl in H or not.
extern int VNL_IN_H; // 25, calculate Vnl in H or not.
extern int ZEEMAN_IN_H; // SunZhiyuan add 


//========================================================================
// EXPLAIN : Parallel information
// GLOBAL VARIABLES :
// NAME : NPROC( global number of process )
// NAME : NPOOL( global number of pools )
// NAME : MY_RANK( global index of process )
// NAME : MY_POOL( global index of pool (count in pool))
// NAME : NPROC_IN_POOL( local number of process in a pool.)
// NAME : RANK_IN_POOL( global index of pool (count in process),
//  	  MY_RANK in each pool)
// NAME : DRANK( index of diag world)
// NAME : DSIZE( number of processors in diag world, only 1 DWORLD exist)
// NAME : DCOLOR( color of each group)
// NAME : GRANK( index of grid world)
// NAME : GSIZE( number of processors in each grid world)
//========================================================================
extern int NPROC;
extern int NPOOL;
extern int MY_RANK;
extern int MY_POOL;
extern int NPROC_IN_POOL;
extern int RANK_IN_POOL;
extern int DRANK;
extern int DSIZE;
extern int DCOLOR;
extern int GRANK;
extern int GSIZE;

//==========================================================
// EXPLAIN : readin file dir, output file ofstream
// GLOBAL VARIABLES :
// NAME : global_in_card
// NAME : global_atom_card
// NAME : global_kpoint_card
// NAME : global_wannier_card
// NAME : global_pseudo_dir
// NAME : global_pseudo_type // mohan add 2013-05-20 (xiaohui add 2013-06-23)
// NAME : global_epm_pseudo_card
// NAME : global_out_dir
// NAME : ofs_running( contain information during runnnig)
// NAME : ofs_warning( contain warning information, including error)
//==========================================================
extern string	global_in_card;
extern string	global_atom_card;
extern string 	global_kpoint_card;
extern string	global_wannier_card;

//#ifndef __EPM
extern string	global_pseudo_dir;
extern string   global_pseudo_type; // mohan add 2013-05-20 (xiaohui add 2013-06-23)
//#else
extern string	global_epm_pseudo_card;
//#endif

extern string 	global_out_dir;

extern ofstream ofs_running;
extern ofstream ofs_warning;

#ifdef __EPM
extern int EPM_SPIN_ORBITAL;
extern bool EPM_ZEEMAN;
extern double EPM_MAG_FIELD_X;
extern double EPM_MAG_FIELD_Y;
extern double EPM_MAG_FIELD_Z;
#endif

//==========================================================
// EXPLAIN : test level for each class
//==========================================================
extern int test_run;

extern int test_input;
extern int test_winput;
extern int test_kpoint;
extern int test_atom;
extern int test_unitcell;
extern int test_symmetry;

extern int test_fft;
extern int test_pw;
extern int test_elec;
extern int test_hm;

extern int test_wf;
extern int test_charge;
extern int test_potential;
extern int test_energy;
extern int test_geo; // mohan add 2011-03-17
//==========================================================
// src_onscaling
//==========================================================
extern int test_atom_arrange;
extern int test_atom_input;
extern int test_grid;
extern int test_grid_driver;
extern int test_overlap;
extern int TEST_FORCE;// mohan add 2011-03-18
extern int test_gridt; // mohan add 2011-03-17
//==========================================================
// src_wannier
//==========================================================
extern int test_spillage;
extern int test_improve_pao;
extern int test_eximport;
extern int test_operation;
extern int test_recon;
extern int test_sph_proj;
extern int test_build;
extern int test_setout;
//==========================================================
// src_epm || src_pseudo
//==========================================================
#ifdef __EPM
extern int test_epm;
extern int test_epm_unitcell;
extern int test_epm_nscf;
#else
extern int test_pseudo_cell;
extern int test_pp;
extern int test_kmesh;
extern int test_mlwf_overlap;
extern int test_mlwf_optimize;
extern int test_ion_dynamics;
#endif
//==========================================================
// src_tools
//==========================================================
extern int test_figure;
extern int test_mathzone;
extern int test_deconstructor;


#endif
