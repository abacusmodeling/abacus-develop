#include "cuda_runtime.h"
#include "cufft.h"
#include "device_launch_parameters.h"

typedef cufftDoubleComplex CUFFT_COMPLEX;

void RoundTrip(const CUFFT_COMPLEX *psi, const double *vr, const int *fft_index, CUFFT_COMPLEX *psic);

