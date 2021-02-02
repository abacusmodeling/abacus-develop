//==========================================================
// AUTHOR : mohan
// DATE : 2021-02-01
//==========================================================
#ifndef RUN_LINEARSCALING_H
#define RUN_LINEARSCALING_H

#include "src_pw/tools.h"
#include "input.h"

class Run_linearscaling
{

public:

    Run_linearscaling();
    ~Run_linearscaling();

	// perform linear-scaling calculations for sparse matrix H and S
	// useless now
    void linear_scaling_line(void);

};

#endif
