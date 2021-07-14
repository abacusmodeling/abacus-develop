//==========================================================
// AUTHOR : liuyu
// DATE : 2021-07-12
//==========================================================
#ifndef RUN_MD_H
#define RUN_MD_H

#include "../src_pw/tools.h"

class Run_md
{

	public:

    Run_md();
    ~Run_md();

    static void md_line(void);
    static void classic_md_line(void);
    static void ai_md_line(void);

};

#endif
