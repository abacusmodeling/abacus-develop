//==========================================================
// Author: Rong Shi
// DATE : 2022-01
//==========================================================

#ifndef DFT_RPA_INTERFACE
#define DFT_RPA_INTERFACE

#include "../input.h"
#include "../module_base/complexmatrix.h"
#include "../module_base/matrix.h"
#include "../src_ri/abfs-vector3_order.h"
#include "src_lcao/local_orbital_charge.h"

#include <string>
#include <vector>

namespace ModuleRPA
{
class DFT_RPA_interface
{
  public:
    DFT_RPA_interface()
    {
        ;
    }
    ~DFT_RPA_interface()
    {
        ;
    }

    void out_for_RPA(Local_Orbital_wfc &lowf);
    void out_eigen_vector(Local_Orbital_wfc &lowf);
    void out_struc();
    void out_bands();
    void out_Cs();
    void out_coulomb_k();
    void print_matrix(char *desc, const ModuleBase::matrix &mat);
    void print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat);
};
} // namespace ModuleRPA
namespace GlobalC
{
extern ModuleRPA::DFT_RPA_interface rpa;
}

#endif
