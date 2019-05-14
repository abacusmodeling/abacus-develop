//=========================================================
//AUTHOR : Peize Lin & jingan
//DATE : 2019-01-16
//=========================================================

#include "wfc_dm_2d.h"
#include "src_global/lapack_connector.h"
#include "src_global/scalapack_connector.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"

#include "src_external/src_test/test_function.h"

void Wfc_Dm_2d::init()
{
	if(GAMMA_ONLY_LOCAL)
	{
		wfc_gamma.resize(NSPIN);
		dm_gamma.resize(NSPIN);
	}
	else
	{
		wfc_k.resize(kv.nks);
		dm_k.resize(kv.nks);
	}
}

void Wfc_Dm_2d::cal_dm()
{
	// dm = wfc.T * wf.wg * wfc.conj()
	if(GAMMA_ONLY_LOCAL)
	{
		for(int is=0; is!=NSPIN; ++is)
		{
			std::vector<double> wg_local(ParaO.ncol,0.0);
			for(int ib_global=0; ib_global!=NBANDS; ++ib_global)
			{
				const int ib_local = ParaO.trace_loc_col[ib_global];
				if(ib_local>=0)
					wg_local[ib_local] = wf.wg(is,ib_global);
			}
			
			// wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
			matrix wg_wfc(wfc_gamma[is]);
			for(int ir=0; ir!=wg_wfc.nr; ++ir)
				LapackConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
			
			// C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
			const double one_float=1.0, zero_float=0.0;
			const int one_int=1;
			const char N_char='N', T_char='T';
			dm_gamma[is].create( wfc_gamma[is].nr, wfc_gamma[is].nc );
			pdgemm_(
				&N_char, &T_char,
				&NLOCAL, &NLOCAL, &NBANDS,
				&one_float,
				wg_wfc.c, &one_int, &one_int, ParaO.desc,
				wfc_gamma[is].c, &one_int, &one_int, ParaO.desc,
				&zero_float,
				dm_gamma[is].c, &one_int, &one_int, ParaO.desc);
		}
	}
	else
	{
		for(int ik=0; ik!=kv.nks; ++ik)
		{			
			std::vector<double> wg_local(ParaO.ncol,0.0);
			for(int ib_global=0; ib_global!=NBANDS; ++ib_global)
			{
				const int ib_local = ParaO.trace_loc_col[ib_global];
				if(ib_local>=0)
					wg_local[ib_local] = wf.wg(ik,ib_global);
			}
			
			// wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw).conj();
			ComplexMatrix wg_wfc = conj(wfc_k[ik]);
			for(int ir=0; ir!=wg_wfc.nr; ++ir)
				LapackConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );

			// C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
			const double one_float=1.0, zero_float=0.0;
			const int one_int=1;
			const char N_char='N', T_char='T';
			dm_k[ik].create( wfc_k[ik].nr, wfc_k[ik].nc );
			pzgemm_(
				&N_char, &T_char,
				&NLOCAL, &NLOCAL, &NBANDS,
				&one_float,
				wg_wfc.c, &one_int, &one_int, ParaO.desc,
				wfc_k[ik].c, &one_int, &one_int, ParaO.desc,
				&zero_float,
				dm_k[ik].c, &one_int, &one_int, ParaO.desc);			
		}
	}
}