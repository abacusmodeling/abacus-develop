//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef RUN_FRAG_H
#define RUN_FRAG_H

#include "src_pw/tools.h"
#include "input.h"

class Run_Frag
{

public:

    Run_Frag();
    ~Run_Frag();

    static void init(void);
	
	// perform Linear Combination of Atomic Orbitals (LCAO) calculations
    void LCAO_line(void);

	// perform plane wave basis calculations
    void plane_wave_line(void);

	// perform linear-scaling calculations for sparse matrix H and S
	// useless now
    void linear_scaling_line(void);

};

#endif
