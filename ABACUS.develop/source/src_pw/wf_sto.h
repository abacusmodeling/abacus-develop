#ifndef WF_STO_H
#define WF_STO_H

#include "tools.h"

class WF_Stochastic 
{
	public:

    WF_Stochastic();
    ~WF_Stochastic();


    ComplexMatrix *evc;  // wavefunctions in the PW basis

    void init_at_1(void);// from init_at_1.f90

	protected:

    //==================================
    // Calculate random wave functions
    // as stochastic wave functions
    //==================================
    void random(ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik)const;

    void check_psi(const ComplexMatrix *psi)const;

};
#endif 
