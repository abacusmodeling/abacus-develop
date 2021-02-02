//==========================================================
// AUTHOR : mohan
// DATE : 2021-02-01
//==========================================================
#ifndef RUN_PW_H
#define RUN_PW_H

#include "src_pw/tools.h"
#include "input.h"

class Run_pw
{

public:

    Run_pw();
    ~Run_pw();

	// perform plane wave basis calculations
    void plane_wave_line(void);

};

#endif
