//==========================================================
// AUTHOR : Peize Lin
// DATE : 2019-06-28
//==========================================================

#ifndef WFC_DM_2D_TEST_H
#define WFC_DM_2D_TEST_H

#include "../../../src_lcao/wfc_dm_2d.h"
#include <iostream>
#include <fstream>

static void os_wfc_2d(ostream &os, const Wfc_Dm_2d & wfc_dm_2d)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		os<<"@@@"<<std::endl;
		for(int is=0; is!=GlobalV::NSPIN; ++is)
		{
			os<<"is:\t"<<is<<std::endl;
			os<<wfc_dm_2d.wfc_gamma[is]<<std::endl;
		}
	}
	else
	{
		os<<"@@@"<<std::endl;
		for(int ik=0; ik!=GlobalC::kv.nks; ++ik)
		{
			os<<"ik:\t"<<ik<<std::endl;
			os<<wfc_dm_2d.dm_k[ik]<<std::endl;
		}
	}
}

static void os_dm_2d(ostream &os, const Wfc_Dm_2d & wfc_dm_2d)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		os<<"@@@"<<std::endl;
		for(int is=0; is!=GlobalV::NSPIN; ++is)
		{
			os<<"is:\t"<<is<<std::endl;
			os<<wfc_dm_2d.dm_gamma[is]<<std::endl;
		}
	}
	else
	{
		os<<"@@@"<<std::endl;
		for(int ik=0; ik!=GlobalC::kv.nks; ++ik)
		{
			os<<"ik:\t"<<ik<<std::endl;
			os<<wfc_dm_2d.dm_k[ik]<<std::endl;
		}
	}
}

static void ofs_wfc_2d(const Wfc_Dm_2d & wfc_dm_2d)
{
	ofstream ofs("wfc_2d_"+TO_STRING(GlobalV::MY_RANK), ofstream::app);
	os_wfc_2d(ofs,wfc_dm_2d);
	ofs.close();
}

static void ofs_dm_2d(const Wfc_Dm_2d & wfc_dm_2d)
{
	ofstream ofs("dm_2d_"+TO_STRING(GlobalV::MY_RANK), ofstream::app);
	os_dm_2d(ofs,wfc_dm_2d);
	ofs.close();
}

#endif