#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "wf_atomic.h"
#include "module_hamilt_general/hamilt.h"

class wavefunc : public WF_atomic
{
	public:

    wavefunc();
    ~wavefunc();

    // allocate memory
    psi::Psi<std::complex<double>>* allocate(const int nks);

    int out_wfc_pw; //qianrui modify 2020-10-19
    int out_wfc_r=0; // Peize Lin add 2021.11.21

    // init_wfc : "random",or "atomic" or "file"
    std::string init_wfc;
	int mem_saver; //1: save evc when doing nscf calculation.
    void wfcinit(psi::Psi<std::complex<double>>* psi_in=nullptr);// from wfcinit.f90
    void wfcinit_k(psi::Psi<std::complex<double>>* psi_in=nullptr);
    int get_starting_nw(void)const;

	// wanf2: save PAO orbitals,
	//void PAO_in_pw_k(const int &ik, ModuleBase::ComplexMatrix &wvf);

	// wanf2: save given localized orbitals. 
	// mohan_to_qianrui: move the LCAO_in_pw_k and LCAO_in_pw_k_q to wavefunc_in_pw.h 
	void LCAO_in_pw_k(const int &ik, ModuleBase::ComplexMatrix &wvf);
	void LCAO_in_pw_k_q(const int &ik, ModuleBase::ComplexMatrix &wvf, ModuleBase::Vector3<double> q);   // pengfei 2016-11-23

	// used if k dependent staff is ready.
	void prepare_k(void);

    int get_R(int ix, int iy, int iz);     // pengfei 2016-11-23

    int iw2it(int iw);
    int iw2ia(int iw);

    void init_after_vc(const int nks); //LiuXh 20180515

    private: // pengfei 2016-11-23

    ModuleBase::Vector3<int> ***R;
    int ** Rmax;
};

namespace hamilt
{

void diago_PAO_in_pw_k2(const int &ik, psi::Psi<std::complex<float>> &wvf, hamilt::Hamilt<float>* phm_in = nullptr);
void diago_PAO_in_pw_k2(const int &ik, psi::Psi<std::complex<double>> &wvf, hamilt::Hamilt<double>* phm_in = nullptr);
void diago_PAO_in_pw_k2(const int &ik, ModuleBase::ComplexMatrix &wvf);

template <typename FPTYPE, typename Device>
void diago_PAO_in_pw_k2(const Device* ctx, const int &ik, psi::Psi<std::complex<FPTYPE>, Device> &wvf, hamilt::Hamilt<FPTYPE, Device>* phm_in = nullptr);

}

#endif //wavefunc
