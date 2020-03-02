//==========================================================

// Author: Lixin He,mohan
// DATE : 2008-11-6
//==========================================================
#ifndef GLOBAL_H
#define GLOBAL_H

#include "src_global/global_variable.h"
#include "src_global/global_function.h"
#include "pw_basis.h"
#include "energy.h"
#include "occupy.h"

#ifndef __EPM

#include "pseudopot_cell_vnl.h"
#include "../run_frag.h"
#include "charge_broyden.h"
#include "potential.h"
#include "electrons.h"
#include "functional.h"

#endif

#include "hamilt.h"
#include "ions.h"
#include "wavefunc.h"
#include "use_fft.h"
#include "klist.h"
//#include "../src_develop/src_wannier/wan_global.h"
#include "output.h"
#include "magnetism.h"
#include "vdwd2.h"
#include "vdwd3.h"	
#include "exx_pw.h"
#include "exx_global.h"
#include "exx_lip.h"

#include "../src_parallel/ft.h"

//==========================================================
// EXPLAIN : define "GLOBAL CLASS"
//==========================================================
extern Use_FFT UFFT;
extern FFT fftwan;
extern kvect kv;
extern output out;

extern PW_Basis pw;
extern energy en;
extern wavefunc wf;
extern Hamilt hm;
extern Exx_pw exxpw;
extern Exx_Global exx_global;
extern Exx_Lip exx_lip;

// if define epm.
#ifdef __EPM

#include "../src_epm/potential_epm.h"
#include "../src_epm/unitcell_epm.h"
#include "../src_epm/so_pw.h"
#include "../src_epm/so_smallbox.h"
#include "../src_epm/zeeman.h"

extern UnitCell_epm ucell;
extern Potential_EPM epm;
extern FFT fftsb;
extern so_pw sopw;
extern so_smallbox sobox;
extern zeeman zm;

// else, plane wave or first principles
#else

#include "symmetry.h"
#include "unitcell_pseudo.h"
#include "../src_parallel/parallel_grid.h"
#include "../src_parallel/parallel_kpoints.h"
extern pseudopot_cell_vnl ppcell;
extern xcfunc xcf;
extern Charge_Broyden chr;
extern potential pot;
extern electrons elec;
//extern wannier wan;
extern Symmetry symm;
extern Magnetism mag;
extern UnitCell_pseudo ucell;
extern Parallel_Grid Pgrid; 
extern Parallel_Kpoints Pkpoints;

extern Vdwd2 vdwd2;		// Peize Lin add 2019-04-26
extern Vdwd3 vdwd3;		// jiyy add 2019-05-18																										 

//zhengdy-soc
#include "soc.h"
extern Soc soc;
#ifdef __FP
#include "../src_lcao/global_fp.h"
#endif

#endif
//==========================================================
//  EXPLAIN : from control_flags.f90
//	EXPLAIN : define "GLOBAL VARIABLES"
//	EXPLAIN : bool variables ,"l" stands for label.
//==========================================================

//extern bool two_fermi_energies;	// TURE : in LSDA calculation case , there are
//extern int isolve;				// "2": Davidson or "1": CG diagonalization
//extern double alpha0;        	// the mixing parameters for the extrapolation
//extern double beta0;         	// ??? of the starting potential
//extern int ngm0;             	// used in mix_rho
//extern int nmix;             	// the number of iteration kept in the history
//extern int iverbosity;       	// type of printing ( 0 few, 1 all )
//extern bool lfixatom;         	// if .TRUE. some atoms is kept fixed
//extern bool lmd;              	// if .TRUE. the calculation is a dynamics
//extern bool lpath;            	// if .TRUE. the calculation is a path optimizations
//extern bool lneb;             	// if .TRUE. the calculation is NEB dynamics
//extern bool lsmd;             	// if .TRUE. the calculation is string dynamics
//extern bool lwf;              	// if .TRUE. the calculation is with wannier functions
//extern bool lphonon;          	// if .TRUE. the calculation is phonon
//extern bool lraman;           	// if .TRUE. the calculation is raman
//extern bool lconstrain;       	// if .TRUE. the calculation is constraint
//extern bool ldamped;          	// if .TRUE. the calculation is a damped dynamics
//extern bool conv_ions;        	// if .TRUE. ionic convergence has been reached
//extern bool noinv;            	// if .TRUE. eliminates inversion symmetry
//extern bool reduce_io;        	// if .TRUE. reduce the I/O to the strict minimum
//extern realArray ns;			//(:,:,:,:), the occupation matrix used in h_psi
//extern realArray nsnew;			//(:,:,:,:) the occupation matrix computed by at
//extern double eth;				// the (corrected) Hubbard contribution
//extern double *Hubbard_U;		//(ntypx),the Hubbard U
//extern int Hubbard_lmax;		// maximum agular momentum of Hubbard states
//extern ComplexMatrix swfcatom;	// = (),
//extern string U_projection;		// 'atomic', 'ortho-atomic', 'file'
//extern double magtot;
//extern double absmag;
//extern int current_spin;
//extern matrix f_inp;
//extern matrix tetra;			// index of k-points in a given tetrahedron
//extern int ntetra;				// specifies how input coordinates are given
//extern bool lberry;
//extern ComplexMatrix becp;		// is used in wf,chr,elec.
//extern bool noncolin;

#endif
