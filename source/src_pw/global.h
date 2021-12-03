#ifndef GLOBAL_H
#define GLOBAL_H

#include "../run_pw.h"
#include "../module_base/global_variable.h"
#include "../module_base/global_function.h"
#include "pw_basis.h"
#include "energy.h"
#include "VNL_in_pw.h"
#include "charge_broyden.h"
#include "potential.h"
#include "xc_type.h"
#include "hamilt.h"
#include "../src_ions/ions.h"
#include "wavefunc.h"
#include "use_fft.h"
#include "klist.h"
#include "../src_io/output.h"
#include "magnetism.h"
#include "vdwd2.h"
#include "vdwd2_parameters.h"
#include "vdwd3.h"
#include "vdwd3_parameters.h"	
#include "../src_io/restart.h" 
#include "exx_global.h"
#include "../src_lcao/exx_lip.h"
#include "../src_parallel/ft.h"

#ifdef __CUDA
#define CHECK_CUDA(func)\
{\
    cudaError_t status = (func);\
    if(status != cudaSuccess)\
    {\
        printf("CUDA API failed at line %d with error: %s (%d)\n",\
            __LINE__, cudaGetErrorString(status), status);\
    }\
}
#endif


//==========================================================
// EXPLAIN : define "GLOBAL CLASS"
//==========================================================
namespace GlobalC
{
extern K_Vectors kv;
extern Use_FFT UFFT;
extern output out;
extern PW_Basis pw;
extern Stochastic_WF sto_wf; //qianrui add 2021-2-5
extern energy en;
extern wavefunc wf;
extern Hamilt hm;
extern Exx_Global exx_global;
extern Exx_Lip exx_lip;
extern pseudopot_cell_vnl ppcell;
}


#include "../module_symmetry/symmetry.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../src_parallel/parallel_grid.h"
#include "../src_parallel/parallel_kpoints.h"
namespace GlobalC
{
extern UnitCell_pseudo ucell;
extern xcfunc xcf;
extern Charge_Broyden CHR;
extern Potential pot;
extern ModuleSymmetry::Symmetry symm;
extern Parallel_Grid Pgrid; 
extern Parallel_Kpoints Pkpoints;
extern Vdwd2_Parameters vdwd2_para;		// Peize Lin add 2021.03.09
extern Vdwd3_Parameters vdwd3_para;		// jiyy add 2021-05-02	
extern Restart restart;	// Peize Lin add 2020.04.04
}

//extern Magnetism mag;

#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif

#endif
