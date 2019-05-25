//=========================================================
//AUTHOR : Peize Lin and jingan
//DATE : 2019-01-16
//=========================================================

#ifndef WFC_DM_2D_H
#define WFC_DM_2D_H

#include "src_global/matrix.h"
#include "src_global/complexmatrix.h"
#include <vector>

class Wfc_Dm_2d
{
public:
	std::vector<matrix> wfc_gamma;			// wfc_gamma[is](ib,iw);
	std::vector<ComplexMatrix> wfc_k;		// wfc_k[ik](ib,iw);
	std::vector<matrix> dm_gamma;			// dm_gamma[is](iw1,iw2);
	std::vector<ComplexMatrix> dm_k;		// dm_k[ik](iw1,iw2);
	
	void init();
	void cal_dm(const matrix &wg);			// wg(ik,ib)
		// dm = wfc.T * wg * wfc.conj()
};

#endif