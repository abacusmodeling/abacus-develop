//==========================================================
// AUTHOR : mohan
// DATE : 2021-01-20
//==========================================================
#ifndef RUN_LCAO_H
#define RUN_LCAO_H

#include "src_pw/tools.h"
#include "input.h"

class Run_lcao
{

	public:

    Run_lcao();
    ~Run_lcao();

	// perform Linear Combination of Atomic Orbitals (LCAO) calculations
    static void lcao_line(void);

};

#endif
