//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#include <vector>
#include <string>
#include <fstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "global_variable.h"
namespace GlobalV
{

//----------------------------------------------------------
// EXPLAIN : Basic Global Variables
// In practice calculation, these values are set in
// input.cpp.
//----------------------------------------------------------
int 	NBANDS = 0;
int		NBANDS_ISTATE = 0; // default number.
int		NLOCAL = 0;	// total number of local basis.

double PSEUDORCUT;
bool PSEUDO_MESH;

std::string	CALCULATION = "scf";
int		EFIELD = 0; // 5: add electric field
int		DIPOLE = 0; // 7: add dipole field

std::string  DFT_FUNCTIONAL = "default";
int 	NSPIN = 1; // LDA
bool	TWO_EFERMI = 0; // two fermi energy, exist only magnetization is fixed.
int 	CURRENT_SPIN = 0;
int 	CURRENT_K = 0;
int		CAL_FORCE = 0;// if cal_force >1, means do the grid integration 'cal_force' times.
double	FORCE_THR = 1.0e-3;
bool	CAL_STRESS = false;
double PRESS1 = 0.0;
double PRESS2 = 0.0;
double PRESS3 = 0.0;
double PRESSURE = 0.0;
std::string	RELAX_METHOD = "bfgs";
std::string  OUT_LEVEL = "ie";
int		RELAX_NMAX = 20;
int     MD_NSTEP = 20;
int 	SCF_NMAX = 50;

std::string	BASIS_TYPE = "pw"; //xiaohui add 2013-09-01
std::string	KS_SOLVER = "cg"; //xiaohui add 2013-09-01
double	SEARCH_RADIUS = -1.0;
bool	SEARCH_PBC = true;
bool	SPARSE_MATRIX = false;

int		DIAGO_PROC = 0;
int 	PW_DIAG_NMAX = 30;
int		DIAGO_CG_PREC = 1; //mohan add 2012-03-31
int 	PW_DIAG_NDIM = 4;
double 	PW_DIAG_THR = 1.0e-2;
int		NB2D = 1;

double 	SCF_THR = 1.0e-9;

std::string	RESTART_MODE = "new";

double DQ = 0.010; // space between Q points of the reciprocal radial tab
int NQX = 10000; // number of points describing reciprocal radial tab

int NURSE = 0; // used for debug.
bool COLOUR = 0;
bool GAMMA_ONLY_LOCAL = 0; // mohan add 2010-10-20
bool GAMMA_ONLY_PW = 0; // mohan add 2012-06-05

int T_IN_H = 1; // mohan add 2010-11-28
int VL_IN_H = 1;
int VNL_IN_H = 1;
int VH_IN_H = 1;
int VION_IN_H = 1;
int ZEEMAN_IN_H = 1;
double  STRESS_THR = 1.0e-2; //LiuXh add 20180515

int ocp=0;
std::string ocp_set = "none";
std::vector<double> ocp_kb(10000);
//int ocp_n=0;
//double ocp_kb[10000];
int  out_mul=0;//qifeng add 2019/9/10
//----------------------------------------------------------
// EXPLAIN : Parallel information
// GLOBAL VARIABLES :
// NAME : NPROC( global number of process )
// NAME : KPAR( global number of pools )
// NAME : MY_RANK( global index of process )
// NAME : MY_POOL( global index of pool (count in pool))
// NAME : NPROC_IN_POOL( local number of process in a pool.)
// NAME : RANK_IN_POOL( global index of pool (count in process),
//  		my_rank in each pool)
//----------------------------------------------------------
int NPROC = 1;
int KPAR = 1;
int MY_RANK = 0;
int MY_POOL = 0;
int NPROC_IN_POOL = 1;
int RANK_IN_POOL = 0;
int DRANK = -1; //mohan add 2012-01-13, must be -1, so we can recognize who didn't in DIAG_WORLD
int DSIZE = KPAR;
int DCOLOR = -1;
int GRANK = MY_RANK;
int GSIZE = DSIZE; 

//----------------------------------------------------------
// EXPLAIN :
//----------------------------------------------------------
std::string	global_in_card = "INPUT";
std::string	stru_file = "STRU";
std::string	global_kpoint_card = "KPT";
std::string	global_wannier_card;

std::string	global_pseudo_dir = "";
std::string global_orbital_dir = "";   // liuyu add 2021-08-14

std::string  global_pseudo_type = "upf"; // mohan add 2013-05-20, default is UPF, we can also use VWR (xiaohui add 2013-06-23)
std::string	global_epm_pseudo_card;
std::string	global_out_dir;
std::string  global_readin_dir; //zhengdy modified

std::ofstream ofs_running;
std::ofstream ofs_warning;

//----------------------------------------------------------
// EXPLAIN : test level for each class
//----------------------------------------------------------
int test_input = 0;
int test_winput = 0;
int test_kpoint = 0;
int test_atom = 0;
int test_unitcell = 0;
int test_symmetry = 0;

int test_fft = 0;
int test_pw = 0;
int test_elec = 0;

int test_wf = 0;
int test_charge = 0;
int test_potential = 0;
int test_energy = 0;
//----------------------------------------------------------
// src_lcao
//----------------------------------------------------------
int test_atom_arrange = 0;
int test_atom_input = 0;
int test_grid = 0;// 4 now
int test_grid_driver = 0;//4 now
int test_overlap = 0;
int TEST_FORCE = 0;//mohan add 2011-03-18
int TEST_STRESS = 0;//zhengdy add 2018-05-16
int test_gridt = 0; // mohan add 2011-03-17
//----------------------------------------------------------
// src_pseudo
//----------------------------------------------------------
int test_pseudo_cell = 0;// 2 : output readin data
int test_pp = 0;// pp: pseudopotential
int test_kmesh = 0;
int test_relax_method = 0;
//----------------------------------------------------------
// src_tools
//----------------------------------------------------------
int test_deconstructor = 0;

//added by zhengdy-soc
bool NONCOLIN = false;
bool LSPINORB = false;
bool DOMAG    = false;
bool DOMAG_Z  = false;
int NPOL      = 1;
int PRENSPIN  = 1;
double soc_lambda = 1.0;

bool FINAL_SCF = false; //LiuXh add 20180619

bool deepks_out_labels = false; //caoyu add 2021-10-16 for DeePKS, wenfei 2022-1-16
bool deepks_scf = false; //caoyu add 2021-10-16 for DeePKS, wenfei 2022-1-16
bool deepks_bandgap = false; //for bandgap label. QO added 2021-12-15
bool deepks_out_unittest = false;

bool deepks_setorb = false;

int vnl_method = 1; //set defauld vnl method as old, added by zhengdy 2021-10-11

bool out_element_info = false; //added by zhengdy 2021-11-26
}
