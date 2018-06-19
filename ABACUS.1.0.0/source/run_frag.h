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

#ifdef __FP
    static void frag_init(void);
    void frag_test(void);
    void frag_pw_line(void);
    void frag_LCAO_line(void);
    void frag_linear_scaling_line(void);
    static void frag_init_after_vc(void); //LiuXh add 20180515
    static void final_calculation_after_vc(void); //LiuXh add 20180619
#endif

    void pw_line(void);


};

#endif
