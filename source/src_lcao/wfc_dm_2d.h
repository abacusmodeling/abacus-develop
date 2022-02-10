//=========================================================
//AUTHOR : Peize Lin and jingan
//DATE : 2019-01-16
//=========================================================

#ifndef WFC_DM_2D_H
#define WFC_DM_2D_H

#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include <vector>

#include "./record_adj.h"   //caoyu add 2022-2-9

class Wfc_Dm_2d
{

	public:

	// wfc stands for wave functions
	std::vector<ModuleBase::matrix> wfc_gamma;			// wfc_gamma[is](ib,iw);
	std::vector<ModuleBase::ComplexMatrix> wfc_k;		// wfc_k[ik](ib,iw);

	// dm stands for density matrix
	std::vector<ModuleBase::matrix> dm_gamma;			// dm_gamma[is](iw1,iw2);
	std::vector<ModuleBase::ComplexMatrix> dm_k;		// dm_k[ik](iw1,iw2);
	
	void init(void);

    // dm = wfc.T * wg * wfc.conj()
    // in multi-k it is dm(k)
	void cal_dm(const ModuleBase::matrix &wg);					// wg(ik,ib), cal all dm 
    // dm(R) = wfc.T * wg * wfc.conj()*kphase, only used in multi-k 
    void cal_dm_R(double** dm2d, Record_adj ra);					// wg(ik,ib), cal dm(R)
};

#endif
