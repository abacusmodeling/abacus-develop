#ifndef TOWannier90_H
#define TOWannier90_H

#include <iostream>
using namespace std;
#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_psi/psi.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#endif

class toWannier90
{
  public:
    // const int k_supercell = 5;
    // const int k_cells = (2 * k_supercell + 1)*(2 * k_supercell + 1)*(2 * k_supercell + 1);
    // const int k_shells = 12;
    // const double large_number = 99999999.0;
    // const double small_number = 0.000001;
    // std::vector<ModuleBase::Vector3<double>> lmn;
    // std::vector<double> dist_shell;
    // std::vector<int> multi;
    // int num_shell_real;
    // int *shell_list_real;
    // double *bweight;

    int num_kpts;
    int cal_num_kpts;
    ModuleBase::Matrix3 recip_lattice;
    std::vector<std::vector<int>> nnlist;
    std::vector<std::vector<ModuleBase::Vector3<double>>> nncell;
    int nntot = 0;
    int num_wannier;
    int *L = nullptr;
    int *m = nullptr;
    int *rvalue = nullptr;
    double *alfa = nullptr;
    ModuleBase::Vector3<double> *R_centre = nullptr;
    std::string wannier_file_name = "seedname";
    int num_exclude_bands = 0;
    int *exclude_bands = nullptr;
    bool *tag_cal_band = nullptr;
    int num_bands;
    bool gamma_only_wannier = false;
    std::string wannier_spin = "up";
    int start_k_index = 0;

    ModuleBase::realArray table_local;
    psi::Psi<std::complex<double>> *unk_inLcao = nullptr;

    toWannier90(int num_kpts, ModuleBase::Matrix3 recip_lattice);
    toWannier90(int num_kpts, ModuleBase::Matrix3 recip_lattice, std::complex<double> ***wfc_k_grid_in);
    ~toWannier90();

    // void kmesh_supercell_sort();
    // void get_nnkpt_first();
    // void kmesh_get_bvectors(int multi, int reference_kpt, double dist_shell,
    // std::vector<ModuleBase::Vector3<double>>& bvector); void get_nnkpt_last();

    void init_wannier(const ModuleBase::matrix& ekb, const psi::Psi<std::complex<double>> *psi = nullptr);
    void read_nnkp();
    void outEIG(const ModuleBase::matrix& ekb);
    void cal_Amn(const psi::Psi<std::complex<double>> &wfc_pw);
    void cal_Mmn(const psi::Psi<std::complex<double>> &wfc_pw);
    void produce_trial_in_pw(const int &ik, ModuleBase::ComplexMatrix &trial_orbitals_k);
    void get_trial_orbitals_lm_k(const int wannier_index,
                                 const int orbital_L,
                                 const int orbital_m,
                                 ModuleBase::matrix &ylm,
                                 ModuleBase::matrix &dr,
                                 ModuleBase::matrix &r,
                                 ModuleBase::matrix &psir,
                                 const int mesh_r,
                                 ModuleBase::Vector3<double> *gk,
                                 const int npw,
                                 ModuleBase::ComplexMatrix &trial_orbitals_k);
    void integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double *table);
    void writeUNK(const psi::Psi<std::complex<double>> &wfc_pw);
    // void ToRealSpace(const int &ik, const int &ib, const ModuleBase::ComplexMatrix *evc, std::complex<double> *psir,
    // const ModuleBase::Vector3<double> G); std::complex<double> unkdotb(const std::complex<double> *psir, const int
    // ikb, const int bandindex, const ModuleBase::ComplexMatrix *wfc_pw);
    std::complex<double> unkdotkb(const int &ik,
                                  const int &ikb,
                                  const int &iband_L,
                                  const int &iband_R,
                                  const ModuleBase::Vector3<double> G,
                                  const psi::Psi<std::complex<double>> &wfc_pw);
    // std::complex<double> gamma_only_cal(const int &ib_L, const int &ib_R, const ModuleBase::ComplexMatrix *wfc_pw,
    // const ModuleBase::Vector3<double> G);

    void lcao2pw_basis(const int ik, ModuleBase::ComplexMatrix &orbital_in_G);
    void getUnkFromLcao();
    void get_lcao_wfc_global_ik(std::complex<double> **ctot, std::complex<double> **cc);

  private:
    std::complex<double> ***wfc_k_grid = nullptr;
};

#endif
