//==========================================================
// AUTHOR : Lixin He, Mohan Chen
// DATE : 2008-11-13
// LAST UPDATE : 2009-03-17 add out_wf
//==========================================================
#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "tools.h"
#include "wf_atomic.h"

class wavefunc : public WF_atomic
{
	public:
    wavefunc();
    ~wavefunc();

    // allocate memory
    void init(const int nks);
#ifdef __FP
    void init_local(void);
#endif

    bool out_wf;

    // et    : (nks,nbnd),eigenvalues of the hamiltonian
    // wg	 : the weight of each k point and band
	double **ekb;		   // band energy at each k point, each band.
	bool   allocate_ekb;   // flag
    matrix wg;

    // start_wfc : "random",or "atomic" or "file"
    string start_wfc;
	int mem_saver; //1: save evc when doing nscf calculation.
    void wfcinit(void);// from wfcinit.f90
    void wfcinit_k();
    int get_starting_nw(void)const;

	// wanf2: save PAO orbitals,
	void PAO_in_pw_k(const int &ik, ComplexMatrix &wvf);

	// wanf2: save given localized orbitals. 
	void LCAO_in_pw_k(const int &ik, ComplexMatrix &wvf);

	// evc: get the initial wave functions from diagnalized the PAO
	// orbitals first.
	void diago_PAO_in_pw_k(const int &ik, ComplexMatrix &wvf);
	// used if k dependent staff is ready.
	void prepare_k(void);
	void diago_PAO_in_pw_k2(const int &ik, ComplexMatrix &wvf);
};

#endif //wavefunc
