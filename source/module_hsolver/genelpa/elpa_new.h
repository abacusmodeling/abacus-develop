#pragma once


extern "C"
{
#include <limits.h>

#include <elpa/elpa_version.h>

struct elpa_struct;
typedef struct elpa_struct *elpa_t;

struct elpa_autotune_struct;
typedef struct elpa_autotune_struct *elpa_autotune_t;

#include <elpa/elpa_constants.h>
#include <elpa/elpa_generated_c_api.h>
// ELPA only provides a C interface header, causing inconsistence of complex
// between C99 (e.g. double complex) and C++11 (std::complex).
// Thus, we have to define a wrapper of complex over the c api
// for compatiability.
#define complex _Complex
#include <elpa/elpa_generated.h>
//#include <elpa/elpa_generic.h>

const char *elpa_strerr(int elpa_error);
}

#define complex _Complex
#include "my_elpa_generic.hpp"
#undef complex