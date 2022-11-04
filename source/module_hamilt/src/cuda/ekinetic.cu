#include "module_hamilt/include/ekinetic.h"
#include <complex>
#include <thrust/complex.h>

using namespace hamilt; 
#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__ void ekinetic_pw(
    const int npw,
    const int max_npw,
    const FPTYPE tpiba2,
    const FPTYPE* gk2,
    thrust::complex<FPTYPE>* hpsi,
    const thrust::complex<FPTYPE>* psi)
{
  const int block_idx = blockIdx.x;
  const int thread_idx = threadIdx.x;
  for (int ii = thread_idx; ii < npw; ii+= blockDim.x) {
    hpsi[block_idx * max_npw + ii] 
      += gk2[ii] * tpiba2 * psi[block_idx * max_npw + ii];
  }
}

template <typename FPTYPE> 
void hamilt::ekinetic_pw_op<FPTYPE, psi::DEVICE_GPU>::operator() (
      const psi::DEVICE_GPU* dev,
      const int& nband,
      const int& npw,
      const int& max_npw,
      const FPTYPE& tpiba2,
      const FPTYPE* gk2_ik,
      std::complex<FPTYPE>* tmhpsi,
      const std::complex<FPTYPE>* tmpsi_in)
{
  // denghui implement 20221019
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ekinetic_pw<FPTYPE><<<nband, THREADS_PER_BLOCK>>>(
    npw, max_npw, tpiba2, // control params
    gk2_ik, // array of data
    reinterpret_cast<thrust::complex<FPTYPE>*>(tmhpsi), // array of data
    reinterpret_cast<const thrust::complex<FPTYPE>*>(tmpsi_in)); // array of data
  // cpu part:
  // for (int ib = 0; ib < nband; ++ib) {
  //   for (int ig = 0; ig < npw; ++ig) {
  //     tmhpsi[ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
  //   }
  //   tmhpsi += max_npw;
  //   tmpsi_in += max_npw;
  // }
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

namespace hamilt{
template struct ekinetic_pw_op<double, psi::DEVICE_GPU>;
}