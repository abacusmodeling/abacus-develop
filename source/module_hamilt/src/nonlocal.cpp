#include "module_hamilt/include/nonlocal.h"

using namespace hamilt;

template <typename FPTYPE> 
struct hamilt::nonlocal_pw_op<FPTYPE, psi::DEVICE_CPU> {
  void operator() (
      const psi::DEVICE_CPU* dev,
      const int& l1,
      const int& l2,
      const int& l3,
      int& sum,
      int& iat,
      const int& spin,
      const int& nkb,
      const int& deeq_x,
      const int& deeq_y,
      const int& deeq_z,
      const FPTYPE* deeq,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp)
  {
    for (int ii = 0; ii < l1; ii++) {
      // each atom has nproj, means this is with structure factor;
      // each projector (each atom) must multiply coefficient
      // with all the other projectors.
      for (int jj = 0; jj < l2; ++jj) 
        for (int kk = 0; kk < l3; kk++) 
          for (int xx = 0; xx < l3; xx++) 
            ps[(sum + ii * l3 + kk) * l2 + jj]
              += deeq[((spin * deeq_x + iat + ii) * deeq_y + xx) * deeq_z + kk] 
              *  becp[jj * nkb + sum + ii * l3 + xx];
    }
    // for (int ww = 0; ww < l1 * l2; ww++) {
    //   const int ii = ww / l2;
    //   const int jj = ww % l2;
    //   for (int tid = 0; tid < l3 * l3; tid++) {
    //     const int kk = tid / l3;
    //     const int xx = tid % l3;
    //     ps[(sum + ii * l3 + kk) * l2 + jj]
    //       += deeq[((spin * deeq_x + iat + ii) * deeq_y + xx) * deeq_z + kk] 
    //       *  becp[jj * nkb + sum + ii * l3 + xx];
    //   }
    // }
    sum += l1 * l3;
    iat += l1;
  }
};

namespace hamilt{
template struct nonlocal_pw_op<double, psi::DEVICE_CPU>;
}