//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-11-6
// LAST MODIFIED : 2021-01-31 by Mohan
//==========================================================

#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------

kvect kv; // mem check in in here.
Use_FFT UFFT; // mohan add 2010-07-22
FFT fftwan;
output out;

PW_Basis pw;
energy en;
wavefunc wf;
Hamilt hm;

UnitCell_pseudo ucell;
pseudopot_cell_vnl ppcell;
xcfunc xcf;
Charge_Broyden CHR;
Magnetism mag;

potential pot;
Symmetry symm;
Parallel_Grid Pgrid; //mohan add 2010-06-06 
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07

Vdwd2 vdwd2(ucell);	// Peize Lin add 2019-04-26
Vdwd3 vdwd3(ucell);	// jiyy add 2019-05-18				

Soc soc; //added by zhengdy-soc														

Restart restart;		// Peize Lin add 2020.04.04

