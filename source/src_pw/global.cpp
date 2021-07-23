#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------
namespace GlobalC
{
K_Vectors kv; // mem check in in here.
Use_FFT UFFT; // mohan add 2010-07-22
output out;
PW_Basis pw;
Stochastic_WF sto_wf;
energy en;
wavefunc wf;
Hamilt hm;
#ifdef __LCAO
Exx_Global exx_global;
Exx_Lip exx_lip(exx_global.info);
#endif
pseudopot_cell_vnl ppcell;
UnitCell_pseudo ucell;
}

xcfunc xcf;
Charge_Broyden CHR;
//Magnetism mag;

Potential pot;
Symmetry symm;
Parallel_Grid Pgrid; //mohan add 2010-06-06 
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07

Vdwd2_Parameters vdwd2_para;	// Peize Lin add 2021-03-09
Vdwd3_Parameters vdwd3_para;	// jiyy add 2021-05-02															

Restart restart; // Peize Lin add 2020.04.04
