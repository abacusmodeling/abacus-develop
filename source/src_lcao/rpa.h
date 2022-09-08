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
#include "local_orbital_charge.h"
//#include "../src_ri/exx_lcao.h"

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

    bool RPA_on=false;

    void out_for_RPA(Local_Orbital_wfc &lowf,Local_Orbital_Charge &loc);
    void out_eigen_vector(Local_Orbital_wfc &lowf);
    void out_struc();
    void out_bands();
    void out_Cs();
    void out_coulomb_k();
    void print_matrix(char *desc, const ModuleBase::matrix &mat);
    void print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat);

    
    void exx_init();
    void exx_cal_ions();
};
} // namespace ModuleRPA
namespace GlobalC
{
extern ModuleRPA::DFT_RPA_interface rpa;

}

#endif
