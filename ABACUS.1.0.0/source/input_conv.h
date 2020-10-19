//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-12-24
//==========================================================
#ifndef INPUT_CONVERT_H
#define INPUT_CONVERT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <vector>

using namespace std;

namespace Input_Conv
{
    void Convert(void);
#ifdef __EPM
    void Convert_EPM(void);
#else
    void Convert_FP(void);
#endif
}


#endif //Input_Convert
