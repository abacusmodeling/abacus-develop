#include "cuda_runtime.h"
#include "cufft.h"
#include "device_launch_parameters.h"

// typedef cufftDoubleComplex CUFFT_COMPLEX;

void RoundTrip_kernel(const float2 *psi, const float *vr, const int *fft_index, float2 *psic);
void RoundTrip_kernel(const double2 *psi, const double *vr, const int *fft_index, double2 *psic);

