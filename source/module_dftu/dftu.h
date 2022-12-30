//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_H
#define DFTU_H

#include "module_cell/unitcell.h"
#include "module_orbital/parallel_orbitals.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_pw/charge_mixing.h"

#include <string>

using namespace std;

//==========================================================
// CLASS :
// NAME : DTFU (DFT+U)
//==========================================================
namespace ModuleDFTU
{

class DFTU
{

  public:
    DFTU(); // constructor
    ~DFTU(); // deconstructor

    //=============================================================
    // In dftu.cpp
    // Initialization & Calculating energy
    //=============================================================
  public:
    // allocate relevant data strcutures
    void init(UnitCell& cell, // unitcell class
              LCAO_Matrix& lm);

    // calculate the energy correction
    void cal_energy_correction(const int istep);
    double get_energy(){return EU;}

    double* U; // U (Hubbard parameter U)
    int* orbital_corr; //
    int omc; // occupation matrix control

  private:
    LCAO_Matrix* LM;
    double EU; //+U energy
    int cal_type = 3; // 1:dftu_tpye=1, dc=1; 2:dftu_type=1, dc=2; 3:dftu_tpye=2, dc=1; 4:dftu_tpye=2, dc=2;
    
    // transform between iwt index and it, ia, L, N and m index
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>
        iatlnmipol2iwt; // iatlnm2iwt[iat][l][n][m][ipol]

    //=============================================================
    // In dftu_hamilt.cpp
    // For calculating contribution to Hamiltonian matrices
    //=============================================================
  public:
    void cal_eff_pot_mat_complex(const int ik, std::complex<double>* eff_pot);
    void cal_eff_pot_mat_real(const int ik, double* eff_pot);
    void cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR);
    void cal_eff_pot_mat_R_complex_double(const int ispin, std::complex<double>* SR, std::complex<double>* HR);

    //=============================================================
    // In dftu_occup.cpp
    // For calculating occupation matrix and saving to locale
    //=============================================================
  public:
    // calculate the local occupation number matrix
    void cal_occup_m_k(const int iter, std::vector<ModuleBase::ComplexMatrix>& dm_k);
    void cal_occup_m_gamma(const int iter, std::vector<ModuleBase::matrix>& dm_gamma);

  private:
    // dftu can be calculated only after locale has been initialed
    bool initialed_locale = false;

    // local occupancy matrix of the correlated subspace
    // locale: the out put local occupation number matrix of correlated electrons in the current electronic step
    // locale_save: the input local occupation number matrix of correlated electrons in the current electronic step
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale; // locale[iat][l][n][spin](m1,m2)
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>
        locale_save; // locale_save[iat][l][n][spin](m1,m2)

    //=============================================================
    // In dftu_tools.cpp
    // For calculating onsite potential, which is used
    // for both Hamiltonian and force/stress
    //=============================================================

    void cal_VU_pot_mat_complex(const int spin, const bool newlocale, std::complex<double>* VU);
    void cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU);

    double get_onebody_eff_pot(const int T,
                               const int iat,
                               const int L,
                               const int N,
                               const int spin,
                               const int m0,
                               const int m1,
                               const bool newlocale);

    //=============================================================
    // In dftu_folding.cpp
    // Subroutines for folding S and dS matrix
    //=============================================================

    void fold_dSR_gamma(const int dim1, const int dim2, double* dSR_gamma);
    // dim = 0 : S, for Hamiltonian
    // dim = 1-3 : dS, for force
    // dim = 4-6 : dS * dR, for stress
    void folding_matrix_k(const int ik, const int dim1, const int dim2, std::complex<double>* mat_k);

    //=============================================================
    // In dftu_force.cpp
    // For calculating force and stress fomr DFT+U
    //=============================================================
  public:
    void force_stress(std::vector<ModuleBase::matrix>& dm_gamma,
                      std::vector<ModuleBase::ComplexMatrix>& dm_k,
                      LCAO_Matrix& lm,
                      ModuleBase::matrix& force_dftu,
                      ModuleBase::matrix& stress_dftu);

  private:
    void cal_force_k(const int ik, const std::complex<double>* rho_VU, ModuleBase::matrix& force_dftu);
    void cal_stress_k(const int ik, const std::complex<double>* rho_VU, ModuleBase::matrix& stress_dftu);
    void cal_force_gamma(const double* rho_VU, ModuleBase::matrix& force_dftu);
    void cal_stress_gamma(const double* rho_VU, ModuleBase::matrix& stress_dftu);

    //=============================================================
    // In dftu_io.cpp
    // For reading/writing/broadcasting/copying relevant data structures
    //=============================================================
  public:
    void output();

  private:
    void write_occup_m(std::ofstream& ofs);
    void read_occup_m(const std::string& fn);
    void local_occup_bcast();
    void copy_locale();

    //=============================================================
    // In dftu_yukawa.cpp
    // Relevant for calculating U using Yukawa potential
    //=============================================================

  public:
    bool Yukawa; // 1:use Yukawa potential; 0: do not use Yukawa potential
    void cal_slater_UJ(double** rho);

  private:
    double lambda; // the parameter in Yukawa potential
    std::vector<std::vector<std::vector<std::vector<double>>>> Fk; // slater integral:Fk[T][L][N][k]
    std::vector<std::vector<std::vector<double>>> U_Yukawa; // U_Yukawa[T][L][N]
    std::vector<std::vector<std::vector<double>>> J_Yukawa; // J_Yukawa[T][L][N]

    void cal_slater_Fk(const int L, const int T); // L:angular momnet, T:atom type
    void cal_yukawa_lambda(double** rho);

    double spherical_Bessel(const int k, const double r, const double lambda);
    double spherical_Hankel(const int k, const double r, const double lambda);
};
} // namespace ModuleDFTU

namespace GlobalC
{
extern ModuleDFTU::DFTU dftu;
}
#endif
