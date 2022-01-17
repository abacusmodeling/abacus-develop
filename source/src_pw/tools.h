#ifndef TOOL_H
#define TOOL_H

using namespace std;
#include <cstdlib>
#include <new>
#include <cassert>

#include <complex>
#include <cmath>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <vector>

#include "../module_base/constants.h"
#include "../module_base/vector3.h"
#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/complexarray.h"
#include "../module_base/complexmatrix.h"

#include "../module_base/memory.h"
#include "../module_base/timer.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/export.h"

#include "../module_base/mathzone.h"
#include "../module_base/mathzone_add1.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/global_file.h"

#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif

#include "../src_parallel/parallel_reduce.h"
#include "../src_parallel/parallel_common.h"

#endif
