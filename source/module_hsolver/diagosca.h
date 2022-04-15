//=====================
// AUTHOR : Peize Lin
// DATE : 2021-11-02
// REFACTORING AUTHOR : Daye Zheng
// DATE : 2022-04-14
//=====================

#ifndef DIAGOSCA_H
#define DIAGOSCA_H

#include "diagh.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include <complex>
#include <utility>
#include <vector>
#include "module_orbital/parallel_orbitals.h"

namespace ModuleHSolver
{

class DiagoSca : public DiagH
{

public:
    void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<double> &psi,
        double *eigenvalue_in)override;

    void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<std::complex<double>> &psi,
        double *eigenvalue_in)override;

private:
	void pdsygvx_diag(
		const int*const desc,
		const int ncol,
		const int nrow,
		const double*const h_mat,
		const double*const s_mat,
		double*const ekb,
		ModulePsi::Psi<double> &wfc_2d);
	void pzhegvx_diag(
		const int*const desc,
		const int ncol,
		const int nrow,
		const std::complex<double>*const h_mat,
		const std::complex<double>*const s_mat,
		double*const ekb,
		ModulePsi::Psi<std::complex<double>> &wfc_2d);

	std::pair<int,std::vector<int>> pdsygvx_once(
		const int*const desc,
		const int ncol,
		const int nrow,
		const double*const h_mat,
		const double*const s_mat,
		double*const ekb,
		ModulePsi::Psi<double> &wfc_2d) const;
	std::pair<int,std::vector<int>> pzhegvx_once(
		const int*const desc,
		const int ncol,
		const int nrow,
		const std::complex<double>*const h_mat,
		const std::complex<double>*const s_mat,
		double*const ekb,
		ModulePsi::Psi<std::complex<double>> &wfc_2d) const;

	int degeneracy_max = 12;			// For reorthogonalized memory. 12 followes siesta.

	void post_processing(const int info, const std::vector<int> &vec);
};

}

#endif