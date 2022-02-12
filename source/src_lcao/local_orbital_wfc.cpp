#include "local_orbital_wfc.h"
#include "../src_pw/global.h"
#include "../src_io/wf_local.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

Local_Orbital_wfc::Local_Orbital_wfc()
{
	allocate_flag = false;
	wfck_flag = false;	
	complex_flag = false;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{

	// used for k-points.
	if(complex_flag && this->wfck_flag)
	{
		for(int i=0; i<GlobalC::kv.nks; i++)
		{
			//for(int j=0; j<GlobalV::NBANDS; j++)
			//{
			//	delete[] this->WFC_K[i][j];
			//}
			delete[] this->WFC_K[i];
			//std::cout<<"delete WFC_K["<<i<<"] success"<<std::endl;
		}
		delete[] this->WFC_K;
		//std::cout<<"delete WFC_K success"<<std::endl;
		if(GlobalV::NLOCAL!= 0 )
		{
			delete[] this->WFC_K_POOL;
			//std::cout<<"delete WFC_K_POOL success"<<std::endl;
		}
	}

}

void Local_Orbital_wfc::allocate_k(const Grid_Technique &gt)
{
	ModuleBase::TITLE("Local_Orbital_wfc","allocate_k");
	if(GlobalV::NLOCAL < GlobalV::NBANDS)
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate","NLOCAL<NBANDS");
	}

	// mohan add the flag 2011-03-02
	// allocate the first part (only once!).
	if(this->wfck_flag == false)
	{
		this->WFC_K = new std::complex<double>**[GlobalC::kv.nks];
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			this->WFC_K[ik] = new std::complex<double>*[GlobalV::NBANDS];
		}
		this->wfck_flag = true;
	}
	
	if(this->complex_flag)
	{
		delete[] this->WFC_K_POOL;
		this->complex_flag = false;
	}
	// allocate the second part.
	//if(gt.lgd != 0) xiaohui modify 2015-02-04, fixed memory bug
	//if(gt.lgd != 0 && this->complex_flag == false)
	if(gt.lgd != 0)
	{
		//std::cout<<"gt.lgd="<<gt.lgd<<" ; GlobalV::NLOCAL="<<GlobalV::NLOCAL<<std::endl; //delete 2015-09-06, xiaohui
		const int page=GlobalV::NBANDS*gt.lgd;
		this->WFC_K_POOL=new std::complex<double> [GlobalC::kv.nks*page];
		ModuleBase::GlobalFunc::ZEROS(WFC_K_POOL, GlobalC::kv.nks*page);
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ib++)
			{
				this->WFC_K[ik][ib] = &WFC_K_POOL[ik*page+ib*gt.lgd];
				//std::cout<<"ik="<<ik<<" ib="<<ib<<std::endl<<"WFC_K address: "<<WFC_K[ik][ib]<<" WFC_K_POOL address: "<<&WFC_K_POOL[ik*page+ib*gt.lgd]<<std::endl;
			}
			//std::cout<<"set WFC_K pointer success, ik: "<<ik<<std::endl;
			ModuleBase::Memory::record("LocalOrbital_Coef","WFC_K",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
			//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MemoryForWaveFunctions (MB)",mem);
			//std::cout<<"WFC_K["<<ik<<"] use "<<mem<<" MB"<<std::endl;
			this->complex_flag = true;
		}
	}

	if(GlobalC::wf.start_wfc == "atomic" )
	{
		
	}
	else if(GlobalC::wf.start_wfc == "file")
	{
		int error;
		std::cout << " Read in wave functions files: " << GlobalC::kv.nkstot << std::endl;
		for(int ik=0; ik<GlobalC::kv.nkstot; ++ik)
		{
			GlobalV::ofs_running << " Read in wave functions " << ik + 1 << std::endl;
			error = WF_Local::read_lowf_complex( this->WFC_K[ik], ik , 0);
		}
#ifdef __MPI
		Parallel_Common::bcast_int(error);
#endif
		GlobalV::ofs_running << " Error=" << error << std::endl;
		if(error==1)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","Can't find the wave function file: LOWF.dat");
		}
		else if(error==2)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, band number doesn't match");
		}
		else if(error==3)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, nlocal doesn't match");
		}
		else if(error==4)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In k-dependent wave function file, k point is not correct");
		}
	}
	else
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc","check the parameter: start_wfc");
	}

	return;
}

