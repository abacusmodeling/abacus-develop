#pragma once

extern "C"
{
#include <elpa/elpa_kernel_constants.h>
#define complex _Complex
#include "my_elpa_generated.h"
#undef complex
}
