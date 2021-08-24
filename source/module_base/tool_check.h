#ifndef TOOL_CHECK_H 
#define TOOL_CHECK_H 

#include <vector>
#include <valarray>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <complex>
#include <cassert>

namespace ModuleBase
{

//==========================================================
// GLOBAL FUNCTIONS :
// Check the input variables
//==========================================================
void CHECK_NAME(std::ifstream &ifs, const std::string &name_in, bool quit=true);
void CHECK_INT(std::ifstream &ifs, const int &v, bool quit=true);
void CHECK_DOUBLE(std::ifstream &ifs, const double &v, bool quit=true);
void CHECK_STRING(std::ifstream &ifs, const std::string &v, bool quit=true);
}

#endif
