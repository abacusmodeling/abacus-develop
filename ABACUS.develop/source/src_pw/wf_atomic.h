#ifndef WF_ATOMIC_H
#define WF_ATOMIC_H

#include "tools.h"
#include "wf_igk.h"

class WF_atomic : public WF_igk
{
	public:
    WF_atomic();
    ~WF_atomic();

	realArray table_local;//mohan add 2009-09-10
    // evc   : (nbnd,npwx),wavefunctions in the PW basis
    // wanf2 : wannier functions in the PW basis
    ComplexMatrix *evc;
    ComplexMatrix *wanf2;

    void init_at_1(void);// from init_at_1.f90
    void print_PAOs(void)const;

	protected:

    void atomic_wfc
    (
        const int ik,
        const int np,
        const int lmax_wfc,
        ComplexMatrix &wfcatom,
        const realArray &table_q,
        const int &table_dimension,
        const double &dq
    )const;

    //==================================
    // Calculate random wave functions
    // as trial wave functions
    //==================================
    void random(ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik)const;
    void check_psi(const ComplexMatrix *psi)const;
};
#endif 
