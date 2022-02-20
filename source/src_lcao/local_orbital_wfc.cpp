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
			//	delete[] this->wfc_k_grid[i][j];
			//}
			delete[] this->wfc_k_grid[i];
			//std::cout<<"delete wfc_k_grid["<<i<<"] success"<<std::endl;
		}
		delete[] this->wfc_k_grid;
		//std::cout<<"delete wfc_k_grid success"<<std::endl;
		if(GlobalV::NLOCAL!= 0 )
		{
			delete[] this->wfc_k_grid2;
			//std::cout<<"delete wfc_k_grid2 success"<<std::endl;
		}
	}

}

void Local_Orbital_wfc::allocate_k(const Grid_Technique& gt,
    std::vector<ModuleBase::ComplexMatrix>& wfc_k)
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
		this->wfc_k_grid = new std::complex<double>**[GlobalC::kv.nks];
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			this->wfc_k_grid[ik] = new std::complex<double>*[GlobalV::NBANDS];
		}
		this->wfck_flag = true;
	}
	
	if(this->complex_flag)
	{
		delete[] this->wfc_k_grid2;
		this->complex_flag = false;
	}
	// allocate the second part.
	//if(gt.lgd != 0) xiaohui modify 2015-02-04, fixed memory bug
	//if(gt.lgd != 0 && this->complex_flag == false)
	if(gt.lgd != 0)
	{
		//std::cout<<"gt.lgd="<<gt.lgd<<" ; GlobalV::NLOCAL="<<GlobalV::NLOCAL<<std::endl; //delete 2015-09-06, xiaohui
		const int page=GlobalV::NBANDS*gt.lgd;
		this->wfc_k_grid2=new std::complex<double> [GlobalC::kv.nks*page];
		ModuleBase::GlobalFunc::ZEROS(wfc_k_grid2, GlobalC::kv.nks*page);
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ib++)
			{
				this->wfc_k_grid[ik][ib] = &wfc_k_grid2[ik*page+ib*gt.lgd];
				//std::cout<<"ik="<<ik<<" ib="<<ib<<std::endl<<"wfc_k_grid address: "<<wfc_k_grid[ik][ib]<<" wfc_k_grid2 address: "<<&wfc_k_grid2[ik*page+ib*gt.lgd]<<std::endl;
			}
			//std::cout<<"set wfc_k_grid pointer success, ik: "<<ik<<std::endl;
			ModuleBase::Memory::record("LocalOrbital_Coef","wfc_k_grid",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
			//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MemoryForWaveFunctions (MB)",mem);
			//std::cout<<"wfc_k_grid["<<ik<<"] use "<<mem<<" MB"<<std::endl;
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
			error = WF_Local::read_lowf_complex( this->wfc_k_grid[ik], ik, &wfc_k);
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

