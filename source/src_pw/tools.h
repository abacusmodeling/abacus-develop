#ifndef TOOL_H
#define TOOL_H

using namespace std;
#include <cstdlib>
#include <new>
#include <cassert>

#include <complex>
#include <cmath>

#include <string>

#include <vector>

#include "../module_base/constants.h"

#include "../module_base/memory.h"
#include "../module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"


#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif

#include "../src_parallel/parallel_reduce.h"
#include "../src_parallel/parallel_common.h"

#endif
