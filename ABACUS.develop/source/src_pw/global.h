//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-11-6
// LAST MODIFIED : 2021-01-31 by Mohan
//==========================================================
#ifndef GLOBAL_H
#define GLOBAL_H

#include "../src_global/global_variable.h"
#include "../src_global/global_function.h"
#include "pw_basis.h"
#include "energy.h"
#include "occupy.h"
#include "pseudopot_cell_vnl.h"
#include "../run_pw.h"
#include "charge_broyden.h"
#include "potential.h"
#include "functional.h"
#include "hamilt.h"
#include "ions.h"
#include "wavefunc.h"
#include "use_fft.h"
#include "klist.h"
#include "../src_pw/output.h"
#include "magnetism.h"
#include "vdwd2.h"
#include "vdwd3.h"				
#include "src_global/restart.h" 
#include "soc.h"

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

#include "symmetry.h"
#include "unitcell_pseudo.h"
#include "../src_parallel/parallel_grid.h"
#include "../src_parallel/parallel_kpoints.h"
extern pseudopot_cell_vnl ppcell;
extern xcfunc xcf;
extern Charge_Broyden CHR;
extern potential pot;
extern Symmetry symm;
extern Magnetism mag;
extern UnitCell_pseudo ucell;
extern Parallel_Grid Pgrid; 
extern Parallel_Kpoints Pkpoints;
extern Vdwd2 vdwd2;		// Peize Lin add 2019-04-26
extern Vdwd3 vdwd3;		// jiyy add 2019-05-18		
extern Restart restart;	// Peize Lin add 2020.04.04  																								 
extern Soc soc; // zhengdy-soc

#include "../src_lcao/global_fp.h"


#endif
